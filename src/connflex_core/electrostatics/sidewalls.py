"""Capacitance between conductor outlines in the layout plane."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from contextlib import AbstractContextManager, nullcontext
from dataclasses import dataclass, field, replace

import numpy as np
from numpy.typing import ArrayLike, NDArray
from shapely import MultiPolygon, Polygon, STRtree

from connflex_core.electrostatics import VACUUM_PERMITTIVITY
from connflex_core.electrostatics.dense import Enter, charge_derivatives, solve_charges
from connflex_core.electrostatics.facing import FacingWalls, facing_plates, feature_sizes
from connflex_core.electrostatics.outlines import Mesh, mesh_outlines
from connflex_core.electrostatics.periodic import (
    Row,
    find_rows,
    interior_cells,
    near_rows,
    solve_rows,
)
from connflex_core.electrostatics.windows import (
    WindowSets,
    cut_cores,
    smallest_window,
    solve_windows,
)
from connflex_core.geometry.heal import as_multipolygon
from connflex_core.params import ElectrostaticParams
from connflex_core.progress import Progress, Stage, stage

Array = NDArray[np.float64]
Index = NDArray[np.intp]
Facing = tuple[Array, Array, Array, FacingWalls]


@dataclass(frozen=True, eq=False)
class SidewallCapacitance:
    """Capacitance between conductor outlines in the layout plane, per meter of thickness.

    Attributes:
        capacitance: n by n matrix in farads per meter: charge per meter of thickness on
            one conductor per volt on another. Symmetric, and every row sums to zero.
        gradient: n by n by m by 3 derivatives of ``capacitance`` for travel of each of m
            rigid bodies along x and y, in farads per square meter, and for its rotation
            about z through its center, in farads per meter and radian.
        hessian: n by n by m by 3 by m by 3 second derivatives of ``capacitance`` for one
            such motion of each of two bodies. The point a body rotates about travels
            with it.
        centers: m by 2 points about which the bodies rotate, in meters.
        feature_size: n by n feature size of each pair of conductors within
            ``interaction_distance``, equation (105): the length of their facing walls
            in units of the gap between them. NaN for other pairs and on the diagonal.
        method: n by n method that found each such pair: ``"dense"``, ``"windowed"``,
            or ``"plates"``; empty for other pairs and on the diagonal. A pair asked to
            be dense whose problem exceeds ``max_elements`` is ``"windowed"``.
        charge: N by n charges on the elements of ``mesh`` that a field solution counts,
            in coulombs per meter of thickness, per volt on each conductor with the
            others at zero; zero elsewhere.
        mesh: Boundary elements on the outlines. Empty when no pair is solved by the
            field method.
        walls: Pieces of wall facing another conductor, which give feature sizes and
            the capacitance of pairs taken as plates.
        windows: Number of field solutions behind the result: 1 for a mesh solved whole,
            the number of windows and rows of identical cells otherwise, counting a
            window that took the solution of a translate and the windows of cells that
            travel took out of rows, 0 when every pair is taken as plates.
        asymmetry: Largest difference between the charge on one conductor per volt on
            another and the same with the two exchanged, over pairs solved by the field
            method, against the largest such entry, before the two were averaged.
            Rounding alone for a mesh solved whole; for windows a measure of what their
            overlap leaves out.
    """

    capacitance: Array
    gradient: Array
    hessian: Array
    centers: Array
    feature_size: Array
    method: NDArray[np.str_]
    charge: Array
    mesh: Mesh
    walls: FacingWalls
    windows: int
    asymmetry: float


@dataclass(frozen=True, eq=False)
class FieldProblem:
    """One field solution: the cores it counts charge on and the pairs it supplies.

    Attributes:
        cores: Elements of each window core.
        pairs: k by 2 pairs of conductors whose capacitance it supplies.
        windows: Elements of each window wherever bodies have taken it, equation (164);
            None before the solution at rest.
        rows: Rows of identical cells whose interior is solved by the periodic kernel
            instead of windows, equation (153).
        loose: Windows of cells that travel has taken out of the interior of their row,
            equation (165).
    """

    cores: list[Index]
    pairs: Index
    windows: WindowSets | None = None
    rows: tuple[Row, ...] = ()
    loose: WindowSets = field(default_factory=WindowSets)


@dataclass(frozen=True, eq=False)
class SidewallPlan:
    """Choices made at rest that the solution keeps when bodies move, equation (106).

    Attributes:
        pairs: k by 2 pairs of conductors within ``interaction_distance``, smaller first.
        feature_size: Feature size of each pair.
        method: Method of each pair after any fallback from dense to windowed.
        mesh: Boundary elements of all outlines; empty when no pair is solved by the
            field method.
        whole: Whether the mesh is solved as one system: every pair is dense and the
            mesh fits ``max_elements``.
        problems: Field solutions in windows, when not ``whole``.
    """

    pairs: Index
    feature_size: Array
    method: tuple[str, ...]
    mesh: Mesh
    whole: bool
    problems: tuple[FieldProblem, ...]

    @property
    def windows(self) -> int:
        """Number of windows and periodic cells the field solutions take."""
        return sum(len(problem.cores) + len(problem.rows) for problem in self.problems)


def sidewall_capacitance(
    conductors: Sequence[Polygon | MultiPolygon],
    params: ElectrostaticParams,
    *,
    bodies: Sequence[int] | None = None,
    centers: ArrayLike | None = None,
    progress: Progress | None = None,
) -> SidewallCapacitance:
    """Capacitance between sidewalls of conductors drawn in the layout plane.

    Each entry of ``conductors`` is one equipotential region, in meters; regions must
    not overlap or touch. Every pair of conductors within ``interaction_distance`` is
    found by the method its feature size selects, equations (104) to (108): below
    ``dense_cutoff`` the field in one dense system, equations (61) to (71), up to
    ``plates_cutoff`` the field in windows, equations (91) to (96), and above it facing
    walls as parallel plates, equations (97) to (103). ``element_size`` is required
    when a pair is solved by the field method.

    ``bodies`` gives the rigid body of each conductor as an index from zero; conductors
    of one body move together. By default every conductor is its own body. ``centers``
    gives the point each body rotates about, by default the centroid of its conductors.
    Gradient and Hessian are the exact derivatives of the computed capacitance for rigid
    motion of the bodies. Stage facing, and when the field is solved stage mesh and then
    assemble, solve, derivatives, and curvature or the one stage windows, log through
    the ``connflex_core`` logger and report to ``progress``.
    """
    regions, body, pivots = _arranged(conductors, bodies, centers)
    with stage("facing", progress) as current:
        facing = _facing(regions, body, pivots, params)
        current.count(pieces=len(facing[3].conductor))
    plan = plan_sidewalls(
        regions, np.arange(len(regions)), facing[3], params, lambda: stage("mesh", progress)
    )
    return solve_sidewalls(
        body,
        pivots,
        plan,
        plan.mesh,
        facing,
        params,
        lambda name, total=None: stage(name, progress, total),
    )[0]


def plan_sidewalls(
    regions: Sequence[MultiPolygon],
    net: Index,
    walls: FacingWalls,
    params: ElectrostaticParams,
    meshing: Callable[[], AbstractContextManager[Stage]] | None = None,
) -> SidewallPlan:
    """Pairs that interact, their feature size and method, and the field problems.

    ``net`` gives the net of every region; pairs of regions take the feature size of
    their pair of nets. ``meshing`` opens the stage in which outlines are meshed when a
    pair is solved by the field method.
    """
    reach = _reach(params)
    first, second = STRtree(list(regions)).query(list(regions), predicate="dwithin", distance=reach)
    keep = first < second
    pairs = np.unique(np.column_stack([first[keep], second[keep]]), axis=0).astype(np.intp)
    sizes = feature_sizes(walls, net)
    feature = np.array(
        [
            sizes.get((min(int(net[a]), int(net[b])), max(int(net[a]), int(net[b]))), 0.0)
            for a, b in pairs
        ]
    )
    method = [
        "dense" if f < params.dense_cutoff else "windowed" if f < params.plates_cutoff else "plates"
        for f in feature
    ]
    field = [k for k, m in enumerate(method) if m != "plates"]
    if not field:
        return SidewallPlan(pairs.reshape(-1, 2), feature, tuple(method), Mesh.empty(), False, ())
    if params.element_size is None:
        a, b = pairs[field[0]]
        raise ValueError(
            f"element_size is required: conductors {a} and {b} have feature size "
            f"{feature[field[0]]:.3g}, below plates_cutoff = {params.plates_cutoff:g}, and are "
            f"solved by the field method"
        )
    size = params.element_size
    with meshing() if meshing is not None else nullcontext() as current:
        mesh = mesh_outlines(regions, size, params.corner_grading)
        if current is not None:
            current.count(elements=len(mesh.conductor))
    every = len(field) == len(pairs) and all(m == "dense" for m in method)
    if every and len(mesh.conductor) <= params.max_elements:
        return SidewallPlan(pairs, feature, tuple(method), mesh, True, ())
    conductors_of: dict[str, set[int]] = {"dense": set(), "windowed": set()}
    for k in field:
        conductors_of[method[k]].update(pairs[k].tolist())
    problems = []
    dense = np.flatnonzero(np.isin(mesh.conductor, sorted(conductors_of["dense"])))
    if dense.size:
        if smallest_window(mesh, dense, size, reach) <= params.max_elements:
            chosen = [k for k in field if method[k] == "dense"]
            problems.append(FieldProblem([dense], pairs[chosen]))
        else:
            method = ["windowed" if m == "dense" else m for m in method]
            conductors_of["windowed"].update(conductors_of["dense"])
    windowed = np.flatnonzero(np.isin(mesh.conductor, sorted(conductors_of["windowed"])))
    if windowed.size:
        chosen = [k for k in field if method[k] == "windowed"]
        rows: list[Row] = []
        if params.periodic_min_cells is not None:
            rows = [
                row
                for row in find_rows(regions, mesh, reach, params.periodic_min_cells)
                if np.isin(mesh.conductor[row.elements], sorted(conductors_of["windowed"])).all()
            ]
            if rows:
                windowed = np.setdiff1d(windowed, np.concatenate([row.elements for row in rows]))
        crowd, groups = None, [windowed] if windowed.size else []
        if rows:

            def crowd(core: Index) -> int:
                return smallest_window(mesh, core, size, reach)

            groups = near_rows(rows, mesh, windowed, reach) if windowed.size else []
        cores = [
            core
            for group in groups
            for core in cut_cores(mesh, params.window_elements, group, crowd)
        ]
        problems.append(FieldProblem(cores, pairs[chosen], rows=tuple(rows)))
    return SidewallPlan(pairs, feature, tuple(method), mesh, False, tuple(problems))


def solve_sidewalls(
    body: Index,
    pivots: Array,
    plan: SidewallPlan,
    mesh: Mesh,
    facing: Facing,
    params: ElectrostaticParams,
    enter: Enter,
    advance: Callable[[], None] | None = None,
) -> tuple[SidewallCapacitance, SidewallPlan]:
    """Capacitance of conductors that move with ``body`` about ``pivots`` by ``plan``.

    ``mesh`` is the mesh of the plan carried with the bodies and ``facing`` the result of
    ``facing_plates`` for the regions where they are. Each pair takes its capacitance
    from its own method, equation (107). ``enter`` opens each stage; with ``advance`` no
    stage is opened for windows and it is called once per window. Also returns the plan
    with the elements of each window, which passed back keep each window as long as it
    holds what its core needs, equation (164).
    """
    count, rigid = len(facing[0]), len(pivots)
    between = np.zeros((count, count))
    slopes = np.zeros((count, count, rigid, 3))
    bends = np.zeros((count, count, rigid, 3, rigid, 3))
    charge = np.zeros((len(mesh.conductor), count))
    uneven, largest, problems, extra = 0.0, 0.0, [], 0
    for k, (a, b) in enumerate(plan.pairs):
        if plan.method[k] == "plates":
            for i, j in ((a, b), (b, a)):
                between[i, j], slopes[i, j], bends[i, j] = (-part[i, j] for part in facing[:3])
    if plan.whole:
        charges = solve_charges(mesh, count, enter)
        solved = [
            (
                charges.capacitance,
                *charge_derivatives(charges, body, pivots, enter),
                charges.charge,
                plan.pairs,
            )
        ]
    else:
        solved = []
        total = plan.windows
        with (
            enter("windows", total=total) if advance is None and total else nullcontext()
        ) as current:
            step = advance or (current.advance if current is not None else (lambda: None))
            for problem in plan.problems:
                *found, kept = solve_windows(
                    mesh,
                    body,
                    pivots,
                    problem.cores,
                    params.element_size or 0.0,
                    params.interaction_distance,
                    params.window_tolerance,
                    params.max_elements,
                    step,
                    problem.windows,
                    problem.pairs,
                )
                if problem.rows:
                    kept_cells = interior_cells(problem.rows, mesh, plan.mesh, body, _reach(params))
                    cores, keys = _loose(problem.rows, kept_cells, params.window_elements)
                    if cores:
                        *part, _ = solve_windows(
                            mesh,
                            body,
                            pivots,
                            cores,
                            params.element_size or 0.0,
                            params.interaction_distance,
                            params.window_tolerance,
                            params.max_elements,
                            lambda: None,
                            problem.loose,
                            problem.pairs,
                            keys,
                        )
                        for summed, piece in zip(found, part, strict=True):
                            summed += piece
                        extra += len(cores)
                    solve_rows(
                        problem.rows,
                        mesh,
                        plan.mesh,
                        body,
                        pivots,
                        _reach(params),
                        (found[0], found[1], found[2], found[3]),
                        kept_cells,
                    )
                    for _ in problem.rows:
                        step()
                solved.append((found[0], found[1], found[2], found[3], problem.pairs))
                problems.append(replace(problem, windows=kept))
    for raw, slope, bend, held, pairs in solved:
        for a, b in pairs:
            for i, j in ((a, b), (b, a)):
                between[i, j] = -(raw[a, b] + raw[b, a]) / 2
                slopes[i, j] = -(slope[a, b] + slope[b, a]) / 2
                bends[i, j] = -(bend[a, b] + bend[b, a]) / 2
            uneven = max(uneven, abs(raw[a, b] - raw[b, a]))
            largest = max(largest, abs(raw[a, b]), abs(raw[b, a]))
        counted = np.any(held != 0, axis=1)
        charge[counted] = held[counted]
    permittivity = VACUUM_PERMITTIVITY * params.ambient_permittivity
    feature = np.full((count, count), np.nan)
    method = np.full((count, count), "", dtype="<U8")
    for k, (a, b) in enumerate(plan.pairs):
        feature[a, b] = feature[b, a] = plan.feature_size[k]
        method[a, b] = method[b, a] = plan.method[k]
    result = SidewallCapacitance(
        capacitance=permittivity * _maxwell(between),
        gradient=permittivity * _maxwell(slopes),
        hessian=permittivity * _maxwell(bends),
        centers=pivots,
        feature_size=feature,
        method=method,
        charge=permittivity * charge,
        mesh=mesh,
        walls=facing[3],
        windows=1 if plan.whole else plan.windows + extra,
        asymmetry=uneven / largest if largest else 0.0,
    )
    return result, replace(plan, problems=tuple(problems)) if problems else plan


def _loose(
    rows: Sequence[Row], kept: Sequence[NDArray[np.bool_]], most: int
) -> tuple[list[Index], list[tuple[int, int, bytes]]]:
    """Cores of the cells that rows no longer count, in blocks of whole cells of at most
    ``most`` elements, each named by its row, its block, and its cells left out."""
    cores, keys = [], []
    for number, (row, mask) in enumerate(zip(rows, kept, strict=True)):
        size = max(1, most // max(1, len(row.cells[0])))
        for block in range(0, len(mask), size):
            out = ~mask[block : block + size]
            if out.any():
                cells = [row.cells[block + int(k)] for k in np.flatnonzero(out)]
                cores.append(np.concatenate(cells))
                keys.append((number, block, out.tobytes()))
    return cores, keys


def _facing(
    regions: Sequence[MultiPolygon], body: Index, pivots: Array, params: ElectrostaticParams
) -> Facing:
    return facing_plates(regions, body, pivots, _reach(params))


def _reach(params: ElectrostaticParams) -> float:
    if params.interaction_distance is None:
        raise ValueError("interaction_distance is required for sidewall capacitance")
    return params.interaction_distance


def _maxwell(between: Array) -> Array:
    """Capacitance matrix from capacitances between pairs, equation (108)."""
    matrix = -between
    diagonal = np.arange(len(between))
    matrix[diagonal, diagonal] = between.sum(axis=1)
    return matrix


def _arranged(
    conductors: Sequence[Polygon | MultiPolygon],
    bodies: Sequence[int] | None,
    centers: ArrayLike | None,
) -> tuple[list[MultiPolygon], NDArray[np.intp], Array]:
    """Regions of conductors, body of each, and point each body rotates about."""
    regions = _disjoint(conductors)
    body = np.arange(len(regions)) if bodies is None else np.asarray(bodies, dtype=np.intp)
    rigid = int(np.max(body, initial=-1)) + 1
    if body.shape != (len(regions),) or not np.array_equal(np.unique(body), np.arange(rigid)):
        raise ValueError("bodies must give every conductor an index, using each index from zero up")
    if centers is None:
        members = [
            [
                polygon
                for region, b in zip(regions, body, strict=True)
                if b == k
                for polygon in region.geoms
            ]
            for k in range(rigid)
        ]
        pivots = np.array([MultiPolygon(group).centroid.coords[0] for group in members])
    else:
        pivots = np.asarray(centers, dtype=float)
        if pivots.shape != (rigid, 2):
            raise ValueError(f"centers must hold one (x, y) point for each of {rigid} bodies")
    return regions, body, pivots


def _disjoint(conductors: Sequence[Polygon | MultiPolygon]) -> list[MultiPolygon]:
    regions = []
    for index, conductor in enumerate(conductors):
        region = as_multipolygon(conductor)
        if region.is_empty or not conductor.is_valid:
            raise ValueError(f"conductor {index} is empty or not a valid polygon")
        regions.append(region)
    polygons = [polygon for region in regions for polygon in region.geoms]
    owner = np.array([index for index, region in enumerate(regions) for _ in region.geoms])
    first, second = STRtree(polygons).query(polygons, predicate="intersects")
    clash = owner[first] < owner[second]
    if clash.any():
        a, b = owner[first][clash][0], owner[second][clash][0]
        raise ValueError(f"conductors {a} and {b} overlap or touch")
    return regions

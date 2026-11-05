"""Capacitance between conductor outlines in the layout plane."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from shapely import MultiPolygon, Polygon, STRtree

from connflex_core.electrostatics import VACUUM_PERMITTIVITY
from connflex_core.electrostatics.dense import Enter, charge_derivatives, solve_charges
from connflex_core.electrostatics.facing import FacingWalls, facing_plates
from connflex_core.electrostatics.outlines import Mesh, mesh_outlines
from connflex_core.electrostatics.windows import cut_cores, solve_windows
from connflex_core.geometry.heal import as_multipolygon
from connflex_core.params import ElectrostaticParams
from connflex_core.progress import Progress, stage

Array = NDArray[np.float64]


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
        charge: N by n charges on the elements of ``mesh``, in coulombs per meter of
            thickness, per volt on each conductor with the others at zero.
        mesh: Boundary elements on the outlines. Method ``"plates"`` has none and
            finds no charge.
        walls: Pieces of wall taken as plates. Method ``"field"`` has none.
        windows: Number of dense solutions behind the result: 1 for a mesh solved whole,
            more for one solved in windows, 0 with method ``"plates"``.
        asymmetry: Largest difference between the charge on one conductor per volt on
            another and the same with the two exchanged, against the largest entry,
            before the two were averaged. Rounding alone for a mesh solved whole; for
            windows a measure of what their overlap leaves out.
    """

    capacitance: Array
    gradient: Array
    hessian: Array
    centers: Array
    charge: Array
    mesh: Mesh
    walls: FacingWalls
    windows: int
    asymmetry: float


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
    not overlap or touch. ``params.sidewall_method`` chooses how. With ``"field"`` the
    two-dimensional field is solved: the set of conductors carries no net charge and
    potential far away floats, equations (61) to (71), and ``element_size`` is required.
    A mesh above ``max_elements`` is solved in windows, equations (91) to (96).
    With ``"plates"`` every wall that faces another conductor within
    ``interaction_distance`` is a parallel plate with it, equations (97) to (103).

    ``bodies`` gives the rigid body of each conductor as an index from zero; conductors
    of one body move together. By default every conductor is its own body. ``centers``
    gives the point each body rotates about, by default the centroid of its conductors.
    Gradient and Hessian are the exact derivatives of the computed capacitance for rigid
    motion of the bodies. Stages mesh, then assemble, solve, derivatives, and curvature
    or the one stage windows, or the one stage facing for plates, log through the
    ``connflex_core`` logger and report to ``progress``.
    """
    if params.sidewall_method is None:
        raise ValueError(
            'sidewall_method is required: "plates" or "field", by the geometry of the device'
        )
    permittivity = VACUUM_PERMITTIVITY * params.ambient_permittivity
    if params.sidewall_method == "plates":
        if params.interaction_distance is None:
            raise ValueError('interaction_distance is required for sidewall_method "plates"')
        regions, body, pivots = _arranged(conductors, bodies, centers)
        with stage("facing", progress) as current:
            matrix, gradient, hessian, walls = facing_plates(
                regions, body, pivots, params.interaction_distance
            )
            current.count(pieces=len(walls.conductor))
        return SidewallCapacitance(
            capacitance=permittivity * matrix,
            gradient=permittivity * gradient,
            hessian=permittivity * hessian,
            centers=pivots,
            charge=np.zeros((0, len(regions))),
            mesh=Mesh.empty(),
            walls=walls,
            windows=0,
            asymmetry=0.0,
        )
    if params.element_size is None:
        raise ValueError('element_size is required for sidewall_method "field"')
    regions, body, pivots = _arranged(conductors, bodies, centers)
    with stage("mesh", progress) as current:
        mesh = mesh_outlines(regions, params.element_size, params.corner_grading)
        current.count(elements=len(mesh.conductor))
    return solve_field(
        mesh, body, pivots, params, lambda name, total=None: stage(name, progress, total)
    )


def solve_field(
    mesh: Mesh,
    body: NDArray[np.intp],
    pivots: Array,
    params: ElectrostaticParams,
    enter: Enter,
    cores: list[NDArray[np.intp]] | None = None,
    advance: Callable[[], None] | None = None,
) -> SidewallCapacitance:
    """Solve the field of a mesh whose conductors move with ``body`` about ``pivots``.

    A mesh of at most ``max_elements`` is one dense solution, equation (70). A larger
    one is solved in windows around ``cores``, by default those of ``cut_cores``,
    equations (91) to (96). ``enter`` opens each stage: ``stage`` bound to a progress
    receiver, or ``quiet`` when the solution is one step of a larger stage. With
    ``advance`` no stage is opened for windows and it is called once per window.
    """
    if params.element_size is None:
        raise ValueError('element_size is required for sidewall_method "field"')
    windows = 1
    if len(mesh.conductor) <= params.max_elements:
        charges = solve_charges(mesh, len(body), enter)
        raw = (
            charges.capacitance,
            *charge_derivatives(charges, body, pivots, enter),
            charges.charge,
        )
    else:
        cores = cut_cores(mesh, params.window_elements) if cores is None else cores
        windows = len(cores)

        def solve(step: Callable[[], None]) -> tuple[Array, Array, Array, Array]:
            return solve_windows(
                mesh,
                body,
                pivots,
                cores,
                params.element_size or 0.0,
                params.interaction_distance,
                params.window_tolerance,
                params.max_elements,
                step,
            )

        if advance is None:
            with enter("windows", total=windows) as current:
                raw = solve(current.advance)
        else:
            raw = solve(advance)
    permittivity = VACUUM_PERMITTIVITY * params.ambient_permittivity
    uneven = np.max(np.abs(raw[0] - raw[0].T)) / np.max(
        np.abs(raw[0]), initial=np.finfo(float).tiny
    )
    return SidewallCapacitance(
        capacitance=permittivity * _paired(raw[0]),
        gradient=permittivity * _paired(raw[1]),
        hessian=permittivity * _paired(raw[2]),
        centers=pivots,
        charge=permittivity * raw[3],
        mesh=mesh,
        walls=FacingWalls.empty(),
        windows=windows,
        asymmetry=float(uneven),
    )


def _paired(raw: Array) -> Array:
    """Symmetric matrix with rows that sum to zero from one whose rows do, equation (94)."""
    mean = (raw + raw.swapaxes(0, 1)) / 2
    diagonal = np.arange(len(raw))
    mean[diagonal, diagonal] = 0.0
    mean[diagonal, diagonal] = -mean.sum(axis=1)
    return mean


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

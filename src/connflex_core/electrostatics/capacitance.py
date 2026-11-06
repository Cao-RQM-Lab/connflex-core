"""Capacitance of a structure at rest."""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import TYPE_CHECKING, Any

import networkx as nx
import numpy as np
import scipy.linalg
import shapely
from numpy.typing import ArrayLike, NDArray
from shapely import MultiPolygon, STRtree, affinity

from connflex_core.electrostatics import VACUUM_PERMITTIVITY
from connflex_core.electrostatics.facing import FacingWalls, facing_plates
from connflex_core.electrostatics.fringe import FringeCapacitance, FringePlan, find_fringe
from connflex_core.electrostatics.outlines import Mesh
from connflex_core.electrostatics.plates import Plate, find_plates
from connflex_core.electrostatics.section import Section
from connflex_core.electrostatics.sidewalls import (
    Facing,
    SidewallCapacitance,
    SidewallPlan,
    plan_sidewalls,
    solve_sidewalls,
)
from connflex_core.geometry.heal import as_multipolygon
from connflex_core.mechanics.assembly import DOF_NAMES, Dof
from connflex_core.params import ElectrostaticParams
from connflex_core.progress import Progress, quiet, stage

if TYPE_CHECKING:
    from connflex_core.model import Structure

Array = NDArray[np.float64]


@dataclass(frozen=True)
class Interaction:
    """How sidewall capacitance between two nets of a cluster is found.

    Attributes:
        nets: The two nets, smaller first; equal for electrodes of one net.
        feature_size: Length of their facing walls in units of the gap between them,
            equation (105). Zero when no wall of one faces the other.
        method: ``"dense"``, ``"windowed"``, or ``"plates"``, equation (106).
    """

    nets: tuple[int, int]
    feature_size: float
    method: str


@dataclass(frozen=True)
class Cluster:
    """Device electrodes whose sidewall capacitance is found together.

    Attributes:
        electrodes: Indices of its electrodes, ascending.
        elements: Number of boundary elements on their outlines, 0 when no pair of its
            nets is solved by the field method.
        windows: Number of field solutions it took: 1 when solved whole, the number of
            windows and rows of identical cells otherwise, counting a window that took
            the solution of a translate and the windows of cells that travel took out
            of rows, 0 when every pair is taken as plates.
        asymmetry: Asymmetry of its field solutions before they were made symmetric,
            as on ``SidewallCapacitance``: for windows a measure of their error.
        walls: Pieces of wall facing another electrode, with electrodes numbered by
            position in ``electrodes``.
        interactions: Feature size and method of each pair of nets whose electrodes in
            the cluster lie within ``interaction_distance`` of each other.
        anchors: Field solutions of it that the surrogate holds, equation (148); 0 when
            ``surrogate_tolerance`` is None or no pair is solved by the field method.
        sections: Distinct cross-sections that the cuts across its walls took for the
            fringe, equation (160); 0 without ``include_fringe``.
    """

    electrodes: tuple[int, ...]
    elements: int
    windows: int
    asymmetry: float
    walls: FacingWalls
    interactions: tuple[Interaction, ...]
    anchors: int = 0
    sections: int = 0


@dataclass(frozen=True, eq=False)
class HessianBlock:
    """Second derivative of capacitance among some electrodes for some degrees of freedom.

    Attributes:
        electrodes: Indices of the k electrodes it couples.
        dofs: Positions in ``CapacitanceModel.dofs`` of its d degrees of freedom.
        values: k by k by d by d second derivatives of the capacitance matrix.
    """

    electrodes: tuple[int, ...]
    dofs: NDArray[np.intp]
    values: Array


@dataclass(frozen=True, eq=False)
class RestGeometry:
    """What a capacitance model keeps as its bodies move.

    Attributes:
        structure: Structure the model was built from.
        params: Parameters it was built with.
        plates: Plates between device and handle electrodes.
        groups: Clusters of device electrodes at rest, and with ``include_fringe`` each
            device electrode alone over the handle that is in no cluster.
    """

    structure: Structure
    params: ElectrostaticParams
    plates: tuple[Plate, ...]
    groups: tuple[RestCluster, ...]


@dataclass(frozen=True, eq=False)
class RestCluster:
    """One cluster of device electrodes at rest.

    Attributes:
        electrodes: Indices of its electrodes, ascending.
        regions: Silicon of each electrode.
        bodies: Bodies its electrodes move with, ascending.
        body: Position in ``bodies`` of the body of each electrode.
        plan: Pairs, feature sizes, methods, mesh, and windows chosen at rest, with the
            sets of elements windows took as bodies traveled, equation (164).
        facing: Result of ``facing_plates`` at rest.
        near: Pairs of its electrodes on different bodies within ``interaction_distance``
            of each other, as positions in ``electrodes``, each with the distance between
            their silicon.
        surrogate: Field solutions kept for reuse at nearby displacements.
        fringe: Handle beneath and cross-sections met, without ``include_fringe`` None.
        beneath: Handle electrodes of ``fringe``.
        held: Bodies of ``fringe``: ``bodies``, then those of handle electrodes that
            move with none of them.
    """

    electrodes: tuple[int, ...]
    regions: tuple[MultiPolygon, ...]
    bodies: tuple[int, ...]
    body: NDArray[np.intp]
    near: tuple[tuple[int, int, float], ...]
    plan: SidewallPlan
    facing: Facing
    surrogate: _Surrogate = field(default_factory=lambda: _Surrogate())
    fringe: FringePlan | None = None
    beneath: tuple[int, ...] = ()
    held: tuple[int, ...] = ()


@dataclass(eq=False)
class _Anchor:
    """One field solution of a cluster, at the motion of its bodies, and how far its
    expansion holds."""

    motion: Array
    found: SidewallCapacitance
    radius: float
    fringe: FringeCapacitance | None = None


@dataclass(eq=False)
class _Surrogate:
    """Field solutions of one cluster at displacements already met, equation (148)."""

    anchors: list[_Anchor] = field(default_factory=list)
    reach: Array | None = None
    start: float = np.inf


@dataclass(frozen=True, eq=False)
class CapacitanceModel:
    """Capacitance between the electrodes of a structure, in farads.

    Every matrix is n by n over ``Structure.electrodes``, symmetric, with rows that sum
    to zero: entry (a, b) is the charge on electrode a per volt on electrode b with all
    others at zero, so ``-matrix[a, b]`` is the capacitance between two electrodes.

    Attributes:
        matrix: Sum of ``oxide``, ``released``, ``sidewall``, and ``fringe``.
        oxide: Part through surviving buried oxide, between device and handle electrodes.
        released: Part across released gaps between device and handle electrodes.
        sidewall: Part between sidewalls of device electrodes in the layout plane.
        fringe: Part of the field above and below the device layer along walls,
            between device electrodes and to the handle beneath, equation (163).
        gradient: n by n by D derivative of ``matrix`` in the degrees of freedom ``dofs``.
        blocks: Second derivative of ``matrix``, one block per plate between two bodies
            and per cluster. Blocks add where they share entries.
        dofs: (body index, degree of freedom) for each of the D degrees of freedom: six
            per movable body in body order, as for the assembled mechanical matrices.
        net: Net of each electrode.
        plates: Plates behind ``oxide`` and ``released``.
        clusters: Groups of two or more device electrodes behind ``sidewall`` and
            ``fringe``.
        unassigned_area: Area in square meters of device silicon in beams and bridges,
            which belongs to no electrode and is not in the model.
        displacement: Displacement of every degree of freedom in ``dofs`` at which the
            model is evaluated, in meters and radians. Zero at rest.
        clearance: Smallest ratio of a gap between two bodies at ``displacement`` to the
            same gap at rest, equation (126): over plates between two bodies, at the
            corner of the plate that comes closest, and over device electrodes of one
            cluster on different bodies, by the distance between their silicon. 1 at
            rest.
        closest: Electrodes on either side of the gap with that ratio. None at rest.
        rest: Plates, clusters, meshes, and windows found at rest, which ``at`` moves.
    """

    matrix: Array
    oxide: Array
    released: Array
    sidewall: Array
    fringe: Array
    gradient: Array
    blocks: tuple[HessianBlock, ...]
    dofs: tuple[Dof, ...]
    net: NDArray[np.intp]
    plates: tuple[Plate, ...]
    clusters: tuple[Cluster, ...]
    unassigned_area: float
    displacement: Array
    clearance: float
    closest: tuple[int, int] | None
    rest: RestGeometry

    def at(
        self, displacement: ArrayLike, *, exact: bool = False, progress: Progress | None = None
    ) -> CapacitanceModel:
        """The model with its bodies displaced, equations (118) and (119).

        ``displacement`` holds one value for each degree of freedom in ``dofs``, measured
        from rest. Every matrix, the gradient, and the second derivative are those at the
        displaced bodies. Clusters, meshes, and cores stay as found at rest; windows and
        rows of identical cells follow travel, equations (164) and (165). With
        ``surrogate_tolerance`` set, clusters solved by the field method come from the
        expansion of an earlier solution where it holds, equations (148) to (150), unless
        ``exact``. Bodies that touch are an error. Stage sidewalls logs through the
        ``connflex_core`` logger and reports to ``progress`` as it does at rest.
        """
        moved = np.asarray(displacement, dtype=float)
        if moved.shape != (len(self.dofs),):
            raise ValueError(
                f"displacement must hold one value for each of {len(self.dofs)} degrees of freedom"
            )
        return _evaluate(self.rest, moved, progress, exact)

    @property
    def net_matrix(self) -> Array:
        """Capacitance matrix between nets, with the electrodes of each net tied together."""
        count = int(np.max(self.net, initial=-1)) + 1
        tied = np.zeros((count, count))
        np.add.at(tied, (self.net[:, None], self.net[None, :]), self.matrix)
        return tied

    def force(self, potentials: ArrayLike) -> Array:
        """Generalized forces on ``dofs`` with ``potentials`` on the electrodes, equation (116)."""
        v = self._volts(potentials)
        force: Array = 0.5 * np.einsum("a,abk,b->k", v, self.gradient, v)
        return force

    def stiffness(self, potentials: ArrayLike) -> Array:
        """D by D electrostatic stiffness with ``potentials`` on the electrodes, equation (116).

        It adds to the mechanical stiffness; a negative diagonal entry softens that motion.
        """
        v = self._volts(potentials)
        stiffness = np.zeros((len(self.dofs), len(self.dofs)))
        for block in self.blocks:
            held = v[list(block.electrodes)]
            stiffness[np.ix_(block.dofs, block.dofs)] -= 0.5 * np.einsum(
                "a,abkl,b->kl", held, block.values, held
            )
        return stiffness

    def to_dict(self) -> dict[str, Any]:
        """JSON-serializable report for review.

        ``pairs`` lists every two electrodes with capacitance between them: its oxide,
        released, sidewall, and fringe parts and their total in farads, their nets, and the
        change of the total with the six degrees of freedom of each body that changes
        it. ``net_matrix``, ``clusters``, ``dofs``, ``displacement``, ``clearance``,
        ``closest``, and ``unassigned_area`` are as on the model.
        """
        pairs = []
        for a, b in np.argwhere(np.triu(self.matrix, 1)).tolist():
            slopes = -self.gradient[a, b].reshape(-1, 6)
            pairs.append(
                {
                    "electrodes": [a, b],
                    "nets": [int(self.net[a]), int(self.net[b])],
                    "oxide": float(-self.oxide[a, b]),
                    "released": float(-self.released[a, b]),
                    "sidewall": float(-self.sidewall[a, b]),
                    "fringe": float(-self.fringe[a, b]),
                    "total": float(-self.matrix[a, b]),
                    "gradient": [
                        {"body": self.dofs[6 * k][0], "values": slope.tolist()}
                        for k, slope in enumerate(slopes)
                        if slope.any()
                    ],
                }
            )
        return {
            "pairs": pairs,
            "net_matrix": self.net_matrix.tolist(),
            "clusters": [
                {
                    "electrodes": list(cluster.electrodes),
                    "elements": cluster.elements,
                    "windows": cluster.windows,
                    "asymmetry": cluster.asymmetry,
                    "interactions": [
                        {
                            "nets": list(interaction.nets),
                            "feature_size": interaction.feature_size,
                            "method": interaction.method,
                        }
                        for interaction in cluster.interactions
                    ],
                }
                for cluster in self.clusters
            ],
            "dofs": [list(dof) for dof in self.dofs],
            "displacement": self.displacement.tolist(),
            "clearance": self.clearance,
            "closest": None if self.closest is None else list(self.closest),
            "unassigned_area": self.unassigned_area,
        }

    def _volts(self, potentials: ArrayLike) -> Array:
        v = np.asarray(potentials, dtype=float)
        if v.shape != (len(self.net),):
            raise ValueError(
                f"potentials must hold one value for each of {len(self.net)} electrodes"
            )
        return v


def find_capacitance(
    structure: Structure,
    params: ElectrostaticParams | None = None,
    *,
    progress: Progress | None = None,
) -> CapacitanceModel:
    """Capacitance between all electrodes of a structure at rest, equations (109) to (116).

    Plates between device and handle electrodes and sidewalls between device electrodes
    are each included unless ``params`` leaves them out. Sidewalls need
    ``interaction_distance``, and ``element_size`` when a pair of nets is solved by the
    field method; each pair of nets in a cluster takes the method its feature size
    selects, equations (104) to (108). With ``include_fringe`` the field above and
    below the device layer along every wall adds to sidewalls, equations (157) to (163),
    for clusters and for device electrodes alone over the handle. Derivatives are in the
    six degrees of freedom of every movable body.
    Stages plates, clusters, and sidewalls log through the
    ``connflex_core`` logger and report to ``progress``; sidewalls advances once per
    cluster, and once per window for a cluster solved in windows.
    """
    params = params or ElectrostaticParams()
    plates = find_plates(structure, params, progress=progress) if params.include_plates else ()
    groups: list[RestCluster] = []
    if params.include_sidewalls:
        reach = params.interaction_distance
        if reach is None:
            raise ValueError(
                "interaction_distance is required for sidewall capacitance; set it or "
                "include_sidewalls=False"
            )
        device = [e for e in structure.electrodes if e.layer == "device"]
        with stage("clusters", progress) as current:
            margin = structure.params.beam_overlap_margin or structure.params.heal_distance
            silicon = [_merged(electrode.region, margin) for electrode in device]
            graph: nx.Graph[int] = nx.Graph()
            graph.add_nodes_from(range(len(device)))
            near = STRtree(silicon).query(silicon, predicate="dwithin", distance=reach)
            graph.add_edges_from(near.T.tolist())
            handle = [e for e in structure.electrodes if e.layer == "handle"]
            below = (
                STRtree([e.region for e in handle]) if handle and params.include_plates else None
            )
            fringe = params.include_fringe
            solved: dict[bytes, Section] = {}

            def kept(group: list[int]) -> bool:
                if len(group) > 1:
                    return True
                return bool(fringe and below is not None and len(below.query(silicon[group[0]])))

            for members in sorted(
                sorted(group) for group in nx.connected_components(graph) if kept(sorted(group))
            ):
                regions = tuple(silicon[k] for k in members)
                bodies = tuple(sorted({device[k].body for k in members}))
                body = np.array([bodies.index(device[k].body) for k in members], np.intp)
                place = {k: position for position, k in enumerate(members)}
                apart = tuple(
                    (place[i], place[j], float(silicon[i].distance(silicon[j])))
                    for i, j in near.T.tolist()
                    if i < j and i in place and j in place and device[i].body != device[j].body
                )
                pivots = np.array([structure.bodies[b].inertial.center[:2] for b in bodies])
                facing = facing_plates(regions, body, pivots, reach)
                nets = np.array([device[k].net for k in members], dtype=np.intp)
                plan, beneath, held = (
                    _fringe_plan(structure, params, regions, bodies, handle, below, solved)
                    if fringe
                    else (None, (), bodies)
                )
                groups.append(
                    RestCluster(
                        electrodes=tuple(device[k].index for k in members),
                        regions=regions,
                        bodies=bodies,
                        body=body,
                        near=apart,
                        plan=plan_sidewalls(regions, nets, facing[3], params),
                        facing=facing,
                        fringe=plan,
                        beneath=beneath,
                        held=held,
                    )
                )
            current.count(
                clusters=sum(len(g.electrodes) > 1 for g in groups),
                windows=sum(g.plan.windows for g in groups),
            )
    rest = RestGeometry(structure, params, plates, tuple(groups))
    return _evaluate(
        rest, np.zeros(6 * sum(not body.anchor for body in structure.bodies)), progress
    )


def _fringe_plan(
    structure: Structure,
    params: ElectrostaticParams,
    regions: tuple[MultiPolygon, ...],
    bodies: tuple[int, ...],
    handle: list[Any],
    below: STRtree | None,
    solved: dict[bytes, Section],
) -> tuple[FringePlan, tuple[int, ...], tuple[int, ...]]:
    """Plan of the fringe of a group: handle electrodes within reach beneath it, their
    bodies among those of the group, and the bodies of the result. Groups share
    ``solved``, the cross-sections met so far."""
    reach = params.interaction_distance or 0.0
    found: list[int] = []
    if below is not None:
        hits = below.query(list(regions), predicate="dwithin", distance=reach)[1]
        found = sorted(set(hits.tolist()))
    held = list(bodies)
    place = []
    for j in found:
        if handle[j].body not in held:
            held.append(handle[j].body)
        place.append(held.index(handle[j].body))
    plan = FringePlan(
        handle=tuple(handle[j].region for j in found),
        handle_body=np.array(place, dtype=np.intp),
        thickness=structure.stack.device.thickness,
        box=structure.stack.box.thickness,
        reach=reach,
        grading=params.corner_grading,
        sections=solved,
    )
    return plan, tuple(handle[j].index for j in found), tuple(held)


def _evaluate(
    rest: RestGeometry, displacement: Array, progress: Progress | None, exact: bool = True
) -> CapacitanceModel:
    """Capacitance model of ``rest`` with its bodies displaced."""
    structure, params = rest.structure, rest.params
    tolerance = params.surrogate_tolerance
    count = len(structure.electrodes)
    oxide, released, sidewall, fringe = (np.zeros((count, count)) for _ in range(4))
    movable = [body.index for body in structure.bodies if not body.anchor]
    first = {index: 6 * position for position, index in enumerate(movable)}
    gradient = np.zeros((count, count, 6 * len(movable)))
    blocks: list[HessianBlock] = []

    def moved(body: int) -> Array:
        return displacement[first[body] : first[body] + 6] if body in first else np.zeros(6)

    clearance, closest = 1.0, None
    for plate in rest.plates:
        value, slope, bend = plate.evaluate(moved(plate.bodies[0]), moved(plate.bodies[1]))
        if displacement.any() and plate.bodies[0] != plate.bodies[1]:
            ratio = plate.lowest(moved(plate.bodies[0]), moved(plate.bodies[1])) / plate.gap
            if ratio < clearance:
                clearance, closest = ratio, plate.electrodes
        pair = np.ix_(plate.electrodes, plate.electrodes)
        oxide[pair] += plate.oxide * _BETWEEN
        released[pair] += (value - plate.oxide) * _BETWEEN
        sides = [side for side, index in enumerate(plate.bodies) if index in first]
        if plate.bodies[0] == plate.bodies[1] or not sides:
            continue
        local = np.concatenate([np.arange(6 * side, 6 * side + 6) for side in sides])
        dofs = np.concatenate([first[plate.bodies[side]] + np.arange(6) for side in sides])
        gradient[np.ix_(plate.electrodes, plate.electrodes, dofs)] += (
            _BETWEEN[:, :, None] * slope[local]
        )
        values = _BETWEEN[:, :, None, None] * bend[np.ix_(local, local)]
        blocks.append(HessianBlock(plate.electrodes, dofs, values))

    clusters: list[Cluster] = []
    groups: list[RestCluster] = []
    thickness = structure.stack.device.thickness
    permittivity = VACUUM_PERMITTIVITY * params.ambient_permittivity
    steps = sum(1 + group.plan.windows for group in rest.groups)
    if rest.groups:
        with stage("sidewalls", progress, total=steps) as current:
            for group in rest.groups:
                rests = np.array([structure.bodies[b].inertial.center for b in group.held])
                levers = [_in_plane(-thickness / 2 - center[2]) for center in rests]
                full = np.array(
                    [
                        np.concatenate([lever @ moved(b), moved(b)[2:5]])
                        for lever, b in zip(levers, group.held, strict=True)
                    ]
                )
                size = len(group.bodies)
                motion = full[:size, :3]
                centers = rests[:size, :2]
                pivots = centers + motion[:, :2]
                regions = _displaced(group, centers, motion)
                for i, j, width in group.near if motion.any() else ():
                    ratio = float(regions[i].distance(regions[j])) / width
                    if ratio < clearance:
                        clearance, closest = ratio, (group.electrodes[i], group.electrodes[j])
                alone = len(group.electrodes) == 1
                fielded = tolerance is not None and any(
                    method != "plates" for method in group.plan.method
                )
                reused = None if exact or not fielded else _reuse(group, full)
                found: SidewallCapacitance | None = None
                shell: FringeCapacitance | None = None
                plan = group.plan
                if reused is not None:
                    found, shell = reused
                    current.advance(group.plan.windows)
                elif not alone:
                    facing = group.facing
                    if motion.any() and "plates" in group.plan.method:
                        facing = facing_plates(
                            regions, group.body, pivots, params.interaction_distance or 0.0
                        )
                    mesh = _carried(group.plan.mesh, group.body, centers, motion)
                    found, plan = solve_sidewalls(
                        group.body, pivots, group.plan, mesh, facing, params, quiet, current.advance
                    )
                if group.fringe is not None and shell is None:
                    shell = find_fringe(
                        group.fringe,
                        regions,
                        group.body,
                        rests[:, :2] + full[:, :2],
                        rests[:, :2],
                        full,
                    )
                if fielded and reused is None and found is not None:
                    _remember(group, full, found, shell, params, structure)
                groups.append(replace(group, plan=plan))
                if found is not None:
                    clusters.append(
                        Cluster(
                            group.electrodes,
                            len(group.plan.mesh.conductor),
                            found.windows,
                            found.asymmetry,
                            found.walls,
                            _interactions(group, structure, plan),
                            len(group.surrogate.anchors),
                            shell.sections if shell is not None else 0,
                        )
                    )
                    held = group.electrodes
                    _add(
                        sidewall,
                        gradient,
                        blocks,
                        held,
                        thickness * found.capacitance,
                        thickness * found.gradient,
                        thickness * found.hessian,
                        [(first.get(b), levers[k]) for k, b in enumerate(group.bodies)],
                    )
                if shell is not None:
                    members = group.electrodes + group.beneath
                    _add(
                        fringe,
                        gradient,
                        blocks,
                        members,
                        permittivity * shell.capacitance,
                        permittivity * shell.gradient,
                        permittivity * shell.hessian,
                        [
                            (first.get(b), np.vstack([levers[k], _RISE]))
                            for k, b in enumerate(group.held)
                        ],
                    )
                current.advance()

    bridges = {link.source for link in structure.links if link.kind == "bridge"}
    return CapacitanceModel(
        matrix=oxide + released + sidewall + fringe,
        oxide=oxide,
        released=released,
        sidewall=sidewall,
        fringe=fringe,
        gradient=gradient,
        blocks=tuple(blocks),
        dofs=tuple((index, name) for index in movable for name in DOF_NAMES),
        net=np.array([electrode.net for electrode in structure.electrodes], dtype=np.intp),
        plates=rest.plates,
        clusters=tuple(clusters),
        unassigned_area=sum(beam.polygon.area for beam in structure.beams)
        + sum(structure.rejected[index].polygon.area for index in bridges),
        displacement=displacement,
        clearance=clearance,
        closest=closest,
        rest=replace(rest, groups=tuple(groups)),
    )


def _add(
    matrix: Array,
    gradient: Array,
    blocks: list[HessianBlock],
    electrodes: tuple[int, ...],
    value: Array,
    slope: Array,
    bend: Array,
    maps: list[tuple[int | None, Array]],
) -> None:
    """Add a part among ``electrodes`` whose derivatives are in local parameters of
    bodies, each mapped to the six degrees of freedom of its body by its matrix; bodies
    without degrees of freedom have None."""
    matrix[np.ix_(electrodes, electrodes)] += value
    free = [k for k, (start, _) in enumerate(maps) if start is not None]
    if not free:
        return
    size, local = len(electrodes), slope.shape[3]
    travel = scipy.linalg.block_diag(*(maps[k][1] for k in free))
    dofs = np.concatenate([int(maps[k][0] or 0) + np.arange(6) for k in free])
    flat = slope[:, :, free].reshape(size, size, local * len(free))
    gradient[np.ix_(electrodes, electrodes, dofs)] += flat @ travel
    curved = bend[:, :, free][:, :, :, :, free].reshape(
        size, size, local * len(free), local * len(free)
    )
    blocks.append(
        HessianBlock(electrodes, dofs, np.einsum("abpq,pk,ql->abkl", curved, travel, travel))
    )


def _reuse(
    group: RestCluster, motion: Array
) -> tuple[SidewallCapacitance, FringeCapacitance | None] | None:
    """Expansion of the nearest solution whose radius reaches ``motion``, equation (149).

    ``motion`` holds six parameters for each body of ``group.held``: travel and turn in
    the plane at device mid-plane, then rise and turns about x and y.
    """
    memo, chosen, best = group.surrogate, None, np.inf
    for anchor in memo.anchors:
        apart = _apart(memo, motion, anchor.motion)
        if apart <= anchor.radius and apart < best:
            chosen, best = anchor, apart
    if chosen is None:
        return None
    return _expand(chosen, motion), _expand_fringe(chosen, motion)


def _remember(
    group: RestCluster,
    motion: Array,
    found: SidewallCapacitance,
    fringe: FringeCapacitance | None,
    params: ElectrostaticParams,
    structure: Structure,
) -> None:
    """Keep a new solution, and set the radius of the nearest kept one from how well
    its expansion predicted this one, equation (150)."""
    memo = group.surrogate
    tolerance = params.surrogate_tolerance or 0.0
    if memo.reach is None:
        centers = [structure.bodies[b].inertial.center[:2] for b in group.bodies]
        reach = np.zeros(len(group.held))
        for k, region in enumerate(group.regions):
            low_x, low_y, high_x, high_y = region.bounds
            corners = np.array([[low_x, low_y], [low_x, high_y], [high_x, low_y], [high_x, high_y]])
            body = int(group.body[k])
            reach[body] = max(reach[body], float(np.max(np.hypot(*(corners - centers[body]).T))))
        memo.reach = reach
        gaps = [width for _, _, width in group.near]
        memo.start = np.sqrt(tolerance) * (
            min(gaps) if gaps else params.interaction_distance or 0.0
        )
    radius = memo.start
    if memo.anchors:
        nearest = min(memo.anchors, key=lambda anchor: _apart(memo, motion, anchor.motion))
        apart = _apart(memo, motion, nearest.motion)
        if 0 < apart <= 4 * nearest.radius:
            guess = _expand(nearest, motion)
            error = max(
                _relative(guess.capacitance, found.capacitance),
                _relative(guess.gradient, found.gradient),
            )
            grown = 2 * nearest.radius
            nearest.radius = (
                min(grown, SAFETY * apart * np.sqrt(tolerance / error)) if error else grown
            )
        radius = nearest.radius
    memo.anchors.append(_Anchor(motion.copy(), found, radius, fringe))


def _apart(memo: _Surrogate, motion: Array, anchor: Array) -> float:
    """Largest travel of a point of any body between two motions, equation (148)."""
    step = motion - anchor
    reach = memo.reach if memo.reach is not None else np.zeros(len(step))
    return float(
        np.max(
            np.hypot(step[:, 0], step[:, 1])
            + abs(step[:, 2]) * reach
            + abs(step[:, 3])
            + (abs(step[:, 4]) + abs(step[:, 5])) * reach
        )
    )


def _expand(anchor: _Anchor, motion: Array) -> SidewallCapacitance:
    """Second-order expansion of a kept solution at ``motion``, equation (149)."""
    found = anchor.found
    size, bodies = found.capacitance.shape[0], found.gradient.shape[2]
    step = (motion[:bodies, :3] - anchor.motion[:bodies, :3]).ravel()
    slope = found.gradient.reshape(size, size, step.size)
    turned = found.hessian.reshape(size, size, step.size, step.size) @ step
    return replace(
        found,
        capacitance=found.capacitance + slope @ step + 0.5 * turned @ step,
        gradient=(slope + turned).reshape(found.gradient.shape),
    )


def _expand_fringe(anchor: _Anchor, motion: Array) -> FringeCapacitance | None:
    """Second-order expansion of the fringe kept with a solution, as equation (149)."""
    found = anchor.fringe
    if found is None:
        return None
    size = found.capacitance.shape[0]
    step = (motion - anchor.motion).ravel()
    slope = found.gradient.reshape(size, size, step.size)
    turned = found.hessian.reshape(size, size, step.size, step.size) @ step
    return FringeCapacitance(
        found.capacitance + slope @ step + 0.5 * turned @ step,
        (slope + turned).reshape(found.gradient.shape),
        found.hessian,
        found.sections,
    )


def _relative(found: Array, truth: Array) -> float:
    scale = float(np.linalg.norm(truth))
    return float(np.linalg.norm(found - truth)) / scale if scale else 0.0


def _turn(angle: float) -> Array:
    """Matrix that turns row vectors by ``angle``."""
    return np.array([[np.cos(angle), np.sin(angle)], [-np.sin(angle), np.cos(angle)]])


def _displaced(group: RestCluster, centers: Array, motion: Array) -> list[MultiPolygon]:
    """Silicon of a cluster with each body moved in the plane, equation (118)."""
    regions = []
    for region, k in zip(group.regions, group.body, strict=True):
        if motion[k].any():
            turn = _turn(motion[k, 2])
            shift = centers[k] + motion[k, :2] - centers[k] @ turn
            matrix = [turn[0, 0], turn[1, 0], turn[0, 1], turn[1, 1], shift[0], shift[1]]
            region = as_multipolygon(affinity.affine_transform(region, matrix))
        regions.append(region)
    first, second = STRtree(regions).query(regions, predicate="intersects")
    touching = group.body[first] != group.body[second]
    if touching.any():
        a, b = first[touching][0], second[touching][0]
        raise ValueError(
            f"electrodes {group.electrodes[a]} and {group.electrodes[b]} touch at this displacement"
        )
    return regions


def _carried(mesh: Mesh, body: NDArray[np.intp], centers: Array, motion: Array) -> Mesh:
    """Mesh of a cluster with the elements of each body moved in the plane."""
    if not motion.any() or not len(mesh.conductor):
        return mesh
    of = body[mesh.conductor]
    ends = []
    for points in (mesh.start, mesh.end):
        carried = np.empty_like(points)
        for k in range(len(centers)):
            mine = of == k
            carried[mine] = (
                centers[k] + motion[k, :2] + (points[mine] - centers[k]) @ _turn(motion[k, 2])
            )
        ends.append(carried)
    return Mesh(ends[0], ends[1], mesh.conductor)


def _interactions(
    group: RestCluster, structure: Structure, plan: SidewallPlan
) -> tuple[Interaction, ...]:
    """Feature size and method of each pair of nets of a cluster."""
    found: dict[tuple[int, int], Interaction] = {}
    for (a, b), size, method in zip(plan.pairs, plan.feature_size, plan.method, strict=True):
        first, second = sorted(structure.electrodes[group.electrodes[k]].net for k in (a, b))
        found.setdefault((first, second), Interaction((first, second), float(size), method))
    return tuple(found[key] for key in sorted(found))


def _merged(region: MultiPolygon, margin: float) -> MultiPolygon:
    """Silicon of an electrode as one region, closed by ``margin`` where its polygons meet."""
    joined = shapely.union_all(shapely.get_parts(region))
    closed = joined.buffer(margin, join_style="mitre").buffer(-margin, join_style="mitre")
    return as_multipolygon(closed)


def _in_plane(lever: float) -> Array:
    """Travel along x and y and rotation about z at device mid-plane, equation (114)."""
    return np.array(
        [
            [1.0, 0.0, 0.0, 0.0, lever, 0.0],
            [0.0, 1.0, 0.0, -lever, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
        ]
    )


_BETWEEN = np.array([[1.0, -1.0], [-1.0, 1.0]])
_RISE = np.eye(6)[2:5]
SAFETY = 0.8

"""Capacitance of a structure at rest."""

from __future__ import annotations

from collections.abc import Callable
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
    interaction_reach,
    plan_sidewalls,
    solve_sidewalls,
)
from connflex_core.electrostatics.translates import taylor
from connflex_core.geometry.heal import as_multipolygon
from connflex_core.mechanics.assembly import DOF_NAMES, Dof
from connflex_core.mechanics.kinematics import point_jacobian
from connflex_core.params import ElectrostaticParams
from connflex_core.progress import Progress, quiet, stage
from connflex_core.workers import Workers

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


@dataclass(eq=False)
class _Beyond:
    """Pairs of device electrodes on different bodies that no cluster measures, with
    their gaps at rest, over every pair within ``distance`` at rest; found once and
    found again at least twice as far out when a motion could bring farther pairs
    within reach."""

    distance: float = 0.0
    pairs: list[tuple[int, int, float]] = field(default_factory=list)


@dataclass(frozen=True, eq=False)
class RestGeometry:
    """What a capacitance model keeps as its bodies move.

    Attributes:
        structure: Structure the model was built from.
        params: Parameters it was built with.
        plates: Plates between device and handle electrodes.
        groups: Clusters of device electrodes at rest, and with ``include_fringe`` each
            device electrode alone over the handle that is in no cluster.
        workers: Threads and processes that run its windows, cross-sections, and
            assembly.
        electrodes: Indices of the device electrodes.
        silicon: Silicon of each of them, as one region.
        radius: Distance from the center of each one's body to the farthest corner of
            its bounding box: how far a turn of the body can carry a point of it.
        paired: Pairs of device electrodes whose gap a cluster measures.
        beyond: Every other pair that a motion so far could have brought within reach,
            with its gap at rest.
    """

    structure: Structure
    params: ElectrostaticParams
    plates: tuple[Plate, ...]
    groups: tuple[RestCluster, ...]
    workers: Workers = field(default_factory=Workers)
    electrodes: tuple[int, ...] = ()
    silicon: tuple[MultiPolygon, ...] = ()
    radius: tuple[float, ...] = ()
    paired: frozenset[tuple[int, int]] = frozenset()
    beyond: _Beyond = field(default_factory=_Beyond)


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
        radius: Distance from the center of each electrode's body to the farthest
            corner of its bounding box: how far a turn of the body can carry a point of
            it.
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
    radius: Array = field(default_factory=lambda: np.zeros(0))
    surrogate: _Surrogate = field(default_factory=lambda: _Surrogate())
    fringe: FringePlan | None = None
    beneath: tuple[int, ...] = ()
    held: tuple[int, ...] = ()


@dataclass(eq=False)
class _Anchor:
    """One field solution of a cluster, at the motion of its bodies, and how far its
    expansion holds; ``block`` and ``fringe_block`` keep its second derivatives mapped
    onto the degrees of freedom, which expansion leaves as they are."""

    motion: Array
    found: SidewallCapacitance
    radius: float
    fringe: FringeCapacitance | None = None
    block: Array | None = None
    fringe_block: Array | None = None


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
        blocks: Second derivative of ``matrix``, among parts with a movable body: one
            block per plate between two bodies, one per cluster for its sidewalls and,
            with ``include_fringe``, one for its fringe, and one per device electrode
            alone over the handle for its fringe. Blocks add where they share entries.
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
            corner of the plate that comes closest, and over device electrodes on
            different bodies within ``interaction_distance`` of each other at rest or at
            ``displacement``, by the distance between their silicon, when that distance
            is set. 1 at rest.
        closest: Electrodes on either side of the gap with that ratio. None at rest and
            wherever no gap is narrower than at rest, when ``clearance`` is 1.
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
        self,
        displacement: ArrayLike,
        *,
        exact: bool = False,
        workers: int | None = None,
        progress: Progress | None = None,
    ) -> CapacitanceModel:
        """The model with its bodies displaced, equations (118) and (119).

        ``displacement`` holds one value for each degree of freedom in ``dofs``, measured
        from rest. Every matrix, the gradient, and the second derivative are those at the
        displaced bodies. Clusters, meshes, and cores stay as found at rest; windows and
        rows of identical cells follow travel, equations (164) and (165). With
        ``surrogate_tolerance`` set, clusters solved by the field method come from the
        expansion of an earlier solution where it holds, equations (148) to (150), unless
        ``exact``. Device electrodes on different bodies that touch, a plate whose gap closes
        at a corner of its overlap, and a displacement at which the model is not finite
        are errors. ``workers`` threads and processes solve the pieces, as many as at
        rest by default.
        Stage sidewalls logs through the ``connflex_core`` logger and reports to
        ``progress`` as it does at rest.
        """
        moved = np.asarray(displacement, dtype=float)
        if moved.shape != (len(self.dofs),):
            raise ValueError(
                f"displacement must hold one value for each of {len(self.dofs)} degrees of freedom"
            )
        model = _evaluate(self.rest, moved, progress, exact, workers)
        finite = np.isfinite(model.matrix).all() and np.isfinite(model.gradient).all()
        if not finite or not all(np.isfinite(block.values).all() for block in model.blocks):
            raise ValueError("capacitance is not finite at this displacement")
        return model

    def expand(self, displacement: ArrayLike) -> CapacitanceModel:
        """The model at a nearby ``displacement`` by expansion of this one, equation
        (149) taken over the whole model.

        ``matrix`` is carried to second order and ``gradient`` to first in the step from
        ``self.displacement``, with the second derivative of ``blocks`` held; the parts
        ``oxide``, ``released``, ``sidewall``, and ``fringe`` stay those of this model.
        ``clearance`` and ``closest`` are those of the geometry at ``displacement``, and
        a plate that closes or electrodes that touch there are errors, as for ``at``.
        Meant for outputs near a displacement the model was evaluated at, where ``at``
        would solve the field again; the error grows with the cube of the step.
        """
        moved = np.asarray(displacement, dtype=float)
        if moved.shape != (len(self.dofs),):
            raise ValueError(
                f"displacement must hold one value for each of {len(self.dofs)} degrees of freedom"
            )
        step = moved - self.displacement
        matrix: Array = self.matrix + self.gradient @ step
        gradient = self.gradient.copy()
        for block in self.blocks:
            dofs = np.asarray(block.dofs)
            local = step[dofs]
            turned: Array = block.values @ local
            gradient[np.ix_(block.electrodes, block.electrodes, dofs)] += turned
            matrix[np.ix_(block.electrodes, block.electrodes)] += 0.5 * turned @ local
        clearance, closest = _clearance(self.rest, moved)
        return replace(
            self,
            matrix=matrix,
            gradient=gradient,
            displacement=moved,
            clearance=clearance,
            closest=closest,
        )

    @property
    def net_matrix(self) -> Array:
        """Capacitance matrix between nets, with the electrodes of each net tied together."""
        count = int(np.max(self.net, initial=-1)) + 1
        tied = np.zeros((count, count))
        np.add.at(tied, (self.net[:, None], self.net[None, :]), self.matrix)
        return tied

    def charge_gradient(self, potentials: ArrayLike) -> Array:
        """Change of the charge on each electrode with each degree of freedom, n by D,
        with ``potentials`` on the electrodes: ``gradient`` contracted with them."""
        v = self._volts(potentials)
        sensitivity: Array = np.einsum("abk,b->ak", self.gradient, v)
        return sensitivity

    def force(self, potentials: ArrayLike) -> Array:
        """Generalized forces on ``dofs`` with ``potentials`` on the electrodes, equation (116)."""
        v = self._volts(potentials)
        force: Array = 0.5 * v @ self.charge_gradient(v)
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
        it. ``clusters`` lists the electrodes, elements, windows, asymmetry, anchors,
        sections, and interactions of each cluster; ``net_matrix``, ``dofs``,
        ``displacement``, ``clearance``, ``closest``, and ``unassigned_area`` are as on
        the model.
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
                    "anchors": cluster.anchors,
                    "sections": cluster.sections,
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
    workers: int | None = None,
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
    ``workers`` threads search windows, cut the fringe, and assemble blocks of rows of
    one dense system, and as many processes solve the field of each window and
    cross-section, every processor by default; sums are taken in a fixed order, so the
    count changes results only within rounding, and the model keeps it for evaluation at
    displaced bodies.
    Stages plates, clusters, and sidewalls log through the
    ``connflex_core`` logger and report to ``progress``; sidewalls advances once per
    cluster, and once per window for a cluster solved in windows.
    """
    params = params or ElectrostaticParams()
    plates = find_plates(structure, params, progress=progress) if params.include_plates else ()
    groups: list[RestCluster] = []
    device = [e for e in structure.electrodes if e.layer == "device"]
    margin = structure.params.beam_overlap_margin
    assert margin is not None  # resolved by ExtractionParams
    silicon = [_merged(electrode.region, margin) for electrode in device]
    paired: set[tuple[int, int]] = set()
    if params.include_sidewalls:
        reach = interaction_reach(params)
        with stage("clusters", progress) as current:
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
                paired.update(
                    (device[members[i]].index, device[members[j]].index) for i, j, _ in apart
                )
                pivots = np.array([structure.bodies[b].inertial.center[:2] for b in bodies])
                radius = np.array(
                    [
                        _corner_reach(region, pivots[b])
                        for region, b in zip(regions, body, strict=True)
                    ]
                )
                facing = facing_plates(regions, body, pivots, reach)
                nets = np.array([device[k].net for k in members], dtype=np.intp)
                narrowest = min((width for _, _, width in apart), default=np.inf)
                plan, beneath, held = (
                    _fringe_plan(
                        structure, params, regions, bodies, handle, below, solved, narrowest
                    )
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
                        radius=radius,
                        fringe=plan,
                        beneath=beneath,
                        held=held,
                    )
                )
            current.count(
                clusters=sum(len(g.electrodes) > 1 for g in groups),
                windows=sum(g.plan.windows for g in groups),
            )
    rest = RestGeometry(
        structure,
        params,
        plates,
        tuple(groups),
        Workers(workers),
        electrodes=tuple(electrode.index for electrode in device),
        silicon=tuple(silicon),
        radius=tuple(
            _corner_reach(region, structure.bodies[electrode.body].inertial.center[:2])
            for region, electrode in zip(silicon, device, strict=True)
        ),
        paired=frozenset(paired),
    )
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
    gap: float,
) -> tuple[FringePlan, tuple[int, ...], tuple[int, ...]]:
    """Plan of the fringe of a group: handle electrodes within reach beneath it, their
    bodies among those of the group, and the bodies of the result. Groups share
    ``solved``, the cross-sections met so far; ``gap`` is the narrowest between
    electrodes of different bodies of the group."""
    reach = interaction_reach(params)
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
        gap=gap,
        sections=solved,
    )
    return plan, tuple(handle[j].index for j in found), tuple(held)


def _evaluate(
    rest: RestGeometry,
    displacement: Array,
    progress: Progress | None,
    exact: bool = True,
    workers: int | None = None,
) -> CapacitanceModel:
    """Capacitance model of ``rest`` with its bodies displaced, by its workers unless
    ``workers`` gives another count."""
    structure, params = rest.structure, rest.params
    pool = rest.workers if workers is None else Workers(workers)
    tolerance = params.surrogate_tolerance
    count = len(structure.electrodes)
    oxide, released, sidewall, fringe = (np.zeros((count, count)) for _ in range(4))
    movable, first, moved = _mover(structure, displacement)
    gradient = np.zeros((count, count, 6 * len(movable)))
    blocks: list[HessianBlock] = []
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
                rests, levers, full = _group_motion(structure, group, moved)
                size = len(group.bodies)
                motion = full[:size, :3]
                centers = rests[:size, :2]
                pivots = centers + motion[:, :2]
                travel = _travel(motion[group.body], group.radius)
                regions = _displaced(group, centers, motion, travel, interaction_reach(params))
                clearance, closest = _nearest(group, regions, travel, clearance, closest)
                alone = len(group.electrodes) == 1
                fielded = tolerance is not None and any(
                    method != "plates" for method in group.plan.method
                )
                reused = None if exact or not fielded else _reuse(group, full)
                found: SidewallCapacitance | None = None
                shell: FringeCapacitance | None = None
                anchor: _Anchor | None = None
                plan = group.plan
                if reused is not None:
                    anchor, found, shell = reused
                    current.advance(group.plan.windows)
                elif not alone:
                    facing = group.facing
                    if motion.any() and "plates" in group.plan.method:
                        facing = facing_plates(
                            regions, group.body, pivots, interaction_reach(params)
                        )
                    mesh = _carried(group.plan.mesh, group.body, centers, motion)
                    found, plan = solve_sidewalls(
                        group.body,
                        pivots,
                        group.plan,
                        mesh,
                        facing,
                        params,
                        quiet,
                        current.advance,
                        pool,
                    )
                if group.fringe is not None and shell is None:
                    shell = find_fringe(
                        group.fringe,
                        regions,
                        group.body,
                        rests[:, :2] + full[:, :2],
                        rests[:, :2],
                        full,
                        pool,
                    )
                if fielded and reused is None and found is not None:
                    anchor = _remember(group, full, found, shell, params, structure)
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
                    kept = anchor.block if anchor is not None else None
                    block = _add(
                        sidewall,
                        gradient,
                        blocks,
                        group.electrodes,
                        thickness * found.capacitance,
                        thickness * found.gradient,
                        None if kept is not None else thickness * found.hessian,
                        [(first.get(b), levers[k]) for k, b in enumerate(group.bodies)],
                        kept,
                    )
                    if anchor is not None:
                        anchor.block = block
                if shell is not None:
                    kept = anchor.fringe_block if anchor is not None else None
                    block = _add(
                        fringe,
                        gradient,
                        blocks,
                        group.electrodes + group.beneath,
                        permittivity * shell.capacitance,
                        permittivity * shell.gradient,
                        None if kept is not None else permittivity * shell.hessian,
                        [
                            (first.get(b), np.vstack([levers[k], _RISE]))
                            for k, b in enumerate(group.held)
                        ],
                        kept,
                    )
                    if anchor is not None:
                        anchor.fringe_block = block
                current.advance()

    if displacement.any() and params.interaction_distance is not None:
        clearance, closest = _beyond_clusters(rest, moved, clearance, closest)
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


def _mover(
    structure: Structure, displacement: Array
) -> tuple[list[int], dict[int, int], Callable[[int], Array]]:
    """Movable bodies, where the six degrees of freedom of each start in
    ``displacement``, and the six of any body, zero for one that is anchored."""
    movable = [body.index for body in structure.bodies if not body.anchor]
    first = {index: 6 * position for position, index in enumerate(movable)}

    def moved(body: int) -> Array:
        return displacement[first[body] : first[body] + 6] if body in first else np.zeros(6)

    return movable, first, moved


def _group_motion(
    structure: Structure, group: RestCluster, moved: Callable[[int], Array]
) -> tuple[Array, list[Array], Array]:
    """Centers at rest of the bodies of a cluster and of its handle, the matrix of each
    that takes its six degrees of freedom to travel and turn in the plane at device
    mid-plane, equation (114), and the six parameters of each body: those three, then
    rise and the turns about x and y."""
    thickness = structure.stack.device.thickness
    rests = np.array([structure.bodies[b].inertial.center for b in group.held])
    levers = [_in_plane(-thickness / 2 - center[2]) for center in rests]
    full = np.array(
        [
            np.concatenate([lever @ moved(b), moved(b)[2:5]])
            for lever, b in zip(levers, group.held, strict=True)
        ]
    )
    return rests, levers, full


def _nearest(
    group: RestCluster,
    regions: list[MultiPolygon],
    travel: Array,
    clearance: float,
    closest: tuple[int, int] | None,
) -> tuple[float, tuple[int, int] | None]:
    """Clearance over the near pairs of a cluster, measured only where ``travel`` leaves
    a gap room to be narrower than ``clearance`` times its width at rest."""
    for i, j, width in group.near:
        if width - travel[i] - travel[j] >= clearance * width:
            continue  # cannot be the narrowest gap so far
        ratio = float(regions[i].distance(regions[j])) / width
        if ratio < clearance:
            clearance, closest = ratio, (group.electrodes[i], group.electrodes[j])
    return clearance, closest


def _clearance(rest: RestGeometry, displacement: Array) -> tuple[float, tuple[int, int] | None]:
    """Clearance of equation (126) at ``displacement`` and the electrodes across its gap,
    from the geometry alone; a plate that closes or electrodes that touch raise as they
    do when the model is evaluated."""
    structure, params = rest.structure, rest.params
    clearance, closest = 1.0, None
    if not displacement.any():
        return clearance, closest
    _, _, moved = _mover(structure, displacement)
    for plate in rest.plates:
        if plate.bodies[0] != plate.bodies[1]:
            lowest = plate.lowest(moved(plate.bodies[0]), moved(plate.bodies[1]))
            if lowest <= 0:
                raise ValueError(
                    f"gap between plates closes: {lowest} m at a corner of the overlap"
                )
            if lowest / plate.gap < clearance:
                clearance, closest = lowest / plate.gap, plate.electrodes
    for group in rest.groups:
        rests, _, full = _group_motion(structure, group, moved)
        size = len(group.bodies)
        motion = full[:size, :3]
        if not motion.any():
            continue
        travel = _travel(motion[group.body], group.radius)
        regions = _displaced(group, rests[:size, :2], motion, travel, interaction_reach(params))
        clearance, closest = _nearest(group, regions, travel, clearance, closest)
    if params.interaction_distance is not None:
        clearance, closest = _beyond_clusters(rest, moved, clearance, closest)
    return clearance, closest


def _add(
    matrix: Array,
    gradient: Array,
    blocks: list[HessianBlock],
    electrodes: tuple[int, ...],
    value: Array,
    slope: Array,
    bend: Array | None,
    maps: list[tuple[int | None, Array]],
    kept: Array | None = None,
) -> Array | None:
    """Add a part among ``electrodes`` whose derivatives are in local parameters of
    bodies, each mapped to the six degrees of freedom of its body by its matrix; bodies
    without degrees of freedom have None. Returns the second derivative mapped onto
    those degrees of freedom, which ``kept`` passes back in place of ``bend`` when the
    part is the expansion of a solution mapped before."""
    matrix[np.ix_(electrodes, electrodes)] += value
    free = [k for k, (start, _) in enumerate(maps) if start is not None]
    if not free:
        return None
    size, local = len(electrodes), slope.shape[3]
    travel = scipy.linalg.block_diag(*(maps[k][1] for k in free))
    dofs = np.concatenate([int(maps[k][0] or 0) + np.arange(6) for k in free])
    flat = slope[:, :, free].reshape(size, size, local * len(free))
    gradient[np.ix_(electrodes, electrodes, dofs)] += flat @ travel
    if kept is None:
        assert bend is not None
        curved = bend[:, :, free][:, :, :, :, free].reshape(
            size, size, local * len(free), local * len(free)
        )
        kept = travel.T @ curved @ travel
    blocks.append(HessianBlock(electrodes, dofs, kept))
    return kept


def _reuse(
    group: RestCluster, motion: Array
) -> tuple[_Anchor, SidewallCapacitance, FringeCapacitance | None] | None:
    """Nearest solution whose radius reaches ``motion`` and its expansion there,
    equation (149).

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
    return chosen, _expand(chosen, motion), _expand_fringe(chosen, motion)


def _remember(
    group: RestCluster,
    motion: Array,
    found: SidewallCapacitance,
    fringe: FringeCapacitance | None,
    params: ElectrostaticParams,
    structure: Structure,
) -> _Anchor:
    """Keep a new solution, and set the radius of the nearest kept one from how well
    its expansion, sidewalls and fringe, predicted this one, equation (150). Returns
    the kept solution that holds the new one.

    A step without motion in the plane can shrink the radius but never grows it: rise
    and tilt leave the sidewalls unchanged and say nothing about their expansion in the
    plane. A solution at a kept motion replaces it.
    """
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
        memo.start = np.sqrt(tolerance) * (min(gaps) if gaps else interaction_reach(params))
    radius = memo.start
    if memo.anchors:
        nearest = min(memo.anchors, key=lambda anchor: _apart(memo, motion, anchor.motion))
        apart = _apart(memo, motion, nearest.motion)
        if apart == 0:  # solved anew where a solution is kept: the new one stands for it
            nearest.found, nearest.fringe = found, fringe
            nearest.block = nearest.fringe_block = None
            return nearest
        if apart <= 4 * nearest.radius:
            guess, shell = _expand(nearest, motion), _expand_fringe(nearest, motion)
            errors = [
                _relative(guess.capacitance, found.capacitance),
                _relative(guess.gradient, found.gradient),
            ]
            if shell is not None and fringe is not None:
                errors += [
                    _relative(shell.capacitance, fringe.capacitance),
                    _relative(shell.gradient, fringe.gradient),
                ]
            error = max(errors)
            bodies = found.gradient.shape[2]
            in_plane = bool((motion[:bodies, :3] != nearest.motion[:bodies, :3]).any())
            grown = 2 * nearest.radius if in_plane else nearest.radius
            nearest.radius = (
                min(grown, SAFETY * apart * np.sqrt(tolerance / error)) if error else grown
            )
        radius = nearest.radius
    anchor = _Anchor(motion.copy(), found, radius, fringe)
    memo.anchors.append(anchor)
    return anchor


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
    value, slope = taylor(
        found.capacitance,
        found.gradient.reshape(size, size, step.size),
        found.hessian.reshape(size, size, step.size, step.size),
        step,
    )
    return replace(found, capacitance=value, gradient=slope.reshape(found.gradient.shape))


def _expand_fringe(anchor: _Anchor, motion: Array) -> FringeCapacitance | None:
    """Second-order expansion of the fringe kept with a solution, as equation (149)."""
    found = anchor.fringe
    if found is None:
        return None
    size = found.capacitance.shape[0]
    step = (motion - anchor.motion).ravel()
    value, slope = taylor(
        found.capacitance,
        found.gradient.reshape(size, size, step.size),
        found.hessian.reshape(size, size, step.size, step.size),
        step,
    )
    return FringeCapacitance(
        value, slope.reshape(found.gradient.shape), found.hessian, found.sections
    )


def _relative(found: Array, truth: Array) -> float:
    scale = float(np.linalg.norm(truth))
    return float(np.linalg.norm(found - truth)) / scale if scale else 0.0


def _turn(angle: float) -> Array:
    """Matrix that turns row vectors by ``angle``."""
    return np.array([[np.cos(angle), np.sin(angle)], [-np.sin(angle), np.cos(angle)]])


def _carry(region: MultiPolygon, center: Array, motion: Array) -> MultiPolygon:
    """``region`` of a body moved in the plane by ``motion`` about ``center``, equation (118)."""
    if not motion.any():
        return region
    turn = _turn(motion[2])
    shift = center + motion[:2] - center @ turn
    matrix = [turn[0, 0], turn[1, 0], turn[0, 1], turn[1, 1], shift[0], shift[1]]
    return as_multipolygon(affinity.affine_transform(region, matrix))


def _travel(motion: Array, radius: Array) -> Array:
    """Farthest any point of each electrode has moved: travel of its body plus its turn
    times ``radius``, for ``motion`` rows of travel and turn of the electrode's body."""
    travel: Array = np.hypot(motion[:, 0], motion[:, 1]) + np.abs(motion[:, 2]) * radius
    return travel


def _corner_reach(region: MultiPolygon, center: Array) -> float:
    """Distance from ``center`` to the farthest corner of the bounding box of ``region``."""
    low_x, low_y, high_x, high_y = region.bounds
    corners = np.array([[low_x, low_y], [low_x, high_y], [high_x, low_y], [high_x, high_y]])
    return float(np.max(np.hypot(*(corners - center).T)))


def _displaced(
    group: RestCluster, centers: Array, motion: Array, travel: Array, reach: float
) -> list[MultiPolygon]:
    """Silicon of a cluster with each body moved in the plane, equation (118), checked
    for electrodes of different bodies that touch unless ``travel``, how far a point of
    each electrode may have moved, rules it out: pairs not in ``near`` were farther than
    ``reach`` apart at rest."""
    regions = [
        _carry(region, centers[k], motion[k])
        for region, k in zip(group.regions, group.body, strict=True)
    ]
    farthest = float(np.max(travel, initial=0.0))
    if 2 * farthest < reach and all(width > travel[i] + travel[j] for i, j, width in group.near):
        return regions
    first, second = STRtree(regions).query(regions, predicate="intersects")
    touching = group.body[first] != group.body[second]
    if touching.any():
        a, b = first[touching][0], second[touching][0]
        raise ValueError(
            f"electrodes {group.electrodes[a]} and {group.electrodes[b]} touch at this displacement"
        )
    return regions


def _beyond_clusters(
    rest: RestGeometry,
    moved: Callable[[int], Array],
    clearance: float,
    closest: tuple[int, int] | None,
) -> tuple[float, tuple[int, int] | None]:
    """Clearance over the gaps that no cluster measures: device electrodes on different
    bodies within reach of each other where the bodies are, equation (126).

    How far a point of each electrode may have moved bounds how far its gaps have
    closed, so only pairs that may have come within reach are carried and measured,
    among those within reach plus twice that bound at rest, which ``rest.beyond`` keeps
    with their gaps at rest.
    """
    structure = rest.structure
    thickness = structure.stack.device.thickness
    bodies = [structure.electrodes[index].body for index in rest.electrodes]
    motions = {}
    centers = []
    for body in bodies:
        center = structure.bodies[body].inertial.center
        centers.append(center[:2])
        if body not in motions:
            motions[body] = _in_plane(-thickness / 2 - center[2]) @ moved(body)
    travel = _travel(np.array([motions[body] for body in bodies]), np.array(rest.radius))
    reach = interaction_reach(rest.params)
    needed = reach + 2 * float(np.max(travel, initial=0.0))
    beyond = rest.beyond
    if needed > beyond.distance:
        needed = max(needed, 2 * beyond.distance)  # so that a motion growing slowly asks once
        first, second = STRtree(rest.silicon).query(
            rest.silicon, predicate="dwithin", distance=needed
        )
        beyond.pairs = [
            (i, j, float(rest.silicon[i].distance(rest.silicon[j])))
            for i, j in zip(first.tolist(), second.tolist(), strict=True)
            if i < j
            and bodies[i] != bodies[j]
            and (rest.electrodes[i], rest.electrodes[j]) not in rest.paired
        ]
        beyond.distance = needed
    carried: dict[int, MultiPolygon] = {}

    def region(k: int) -> MultiPolygon:
        if k not in carried:
            carried[k] = _carry(rest.silicon[k], centers[k], motions[bodies[k]])
        return carried[k]

    for i, j, at_rest in beyond.pairs:
        if at_rest - travel[i] - travel[j] > reach:
            continue  # cannot have come within reach
        gap = float(region(i).distance(region(j)))
        pair = (rest.electrodes[i], rest.electrodes[j])
        if gap <= 0:
            raise ValueError(f"electrodes {pair[0]} and {pair[1]} touch at this displacement")
        ratio = gap / at_rest
        if gap <= reach and ratio < clearance:
            clearance, closest = ratio, pair
    return clearance, closest


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
    point = point_jacobian([0.0, 0.0, lever])
    return np.vstack([point[:2], [[0.0, 0.0, 0.0, 0.0, 0.0, 1.0]]])


_BETWEEN = np.array([[1.0, -1.0], [-1.0, 1.0]])
_RISE = np.eye(6)[2:5]
SAFETY = 0.5

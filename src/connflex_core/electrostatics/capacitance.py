"""Capacitance of a structure at rest."""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import TYPE_CHECKING, Any

import networkx as nx
import numpy as np
import scipy.linalg
import shapely
from numpy.typing import ArrayLike, NDArray
from shapely import MultiPolygon, STRtree, affinity

from connflex_core.electrostatics import VACUUM_PERMITTIVITY
from connflex_core.electrostatics.facing import FacingWalls, facing_plates
from connflex_core.electrostatics.outlines import Mesh, mesh_outlines
from connflex_core.electrostatics.plates import Plate, find_plates
from connflex_core.electrostatics.sidewalls import solve_field
from connflex_core.electrostatics.windows import cut_cores
from connflex_core.geometry.heal import as_multipolygon
from connflex_core.mechanics.assembly import DOF_NAMES, Dof
from connflex_core.params import ElectrostaticParams
from connflex_core.progress import Progress, quiet, stage

if TYPE_CHECKING:
    from connflex_core.model import Structure

Array = NDArray[np.float64]


@dataclass(frozen=True)
class Cluster:
    """Device electrodes whose sidewall capacitance is found together.

    Attributes:
        electrodes: Indices of its electrodes, ascending.
        elements: Number of boundary elements on their outlines with method ``"field"``,
            of pieces of wall taken as plates with method ``"plates"``.
        windows: Number of dense solutions it took: 1 when solved whole, more when
            solved in windows, 0 with method ``"plates"``.
        asymmetry: Asymmetry of its capacitance matrix before it was made symmetric,
            as on ``SidewallCapacitance``: for windows a measure of their error.
        walls: Pieces of wall taken as plates, with conductors numbered by position in
            ``electrodes``. Method ``"field"`` has none.
    """

    electrodes: tuple[int, ...]
    elements: int
    windows: int
    asymmetry: float
    walls: FacingWalls


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
        groups: Clusters of device electrodes at rest.
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
        mesh: Boundary elements on the outlines with method ``"field"``, else None.
        cores: Cores of its windows; empty when it is one dense solution.
        windows: Elements of each window as chosen at rest; None before that.
    """

    electrodes: tuple[int, ...]
    regions: tuple[MultiPolygon, ...]
    bodies: tuple[int, ...]
    body: NDArray[np.intp]
    mesh: Mesh | None
    cores: list[NDArray[np.intp]]
    windows: list[NDArray[np.intp]] | None


@dataclass(frozen=True, eq=False)
class CapacitanceModel:
    """Capacitance between the electrodes of a structure, in farads.

    Every matrix is n by n over ``Structure.electrodes``, symmetric, with rows that sum
    to zero: entry (a, b) is the charge on electrode a per volt on electrode b with all
    others at zero, so ``-matrix[a, b]`` is the capacitance between two electrodes.

    Attributes:
        matrix: Sum of ``oxide``, ``released``, and ``sidewall``.
        oxide: Part through surviving buried oxide, between device and handle electrodes.
        released: Part across released gaps between device and handle electrodes.
        sidewall: Part between sidewalls of device electrodes.
        gradient: n by n by D derivative of ``matrix`` in the degrees of freedom ``dofs``.
        blocks: Second derivative of ``matrix``, one block per plate between two bodies
            and per cluster. Blocks add where they share entries.
        dofs: (body index, degree of freedom) for each of the D degrees of freedom: six
            per movable body in body order, as for the assembled mechanical matrices.
        net: Net of each electrode.
        plates: Plates behind ``oxide`` and ``released``.
        clusters: Groups of two or more device electrodes behind ``sidewall``.
        unassigned_area: Area in square meters of device silicon in beams and bridges,
            which belongs to no electrode and is not in the model.
        displacement: Displacement of every degree of freedom in ``dofs`` at which the
            model is evaluated, in meters and radians. Zero at rest.
        rest: Plates, clusters, meshes, and windows found at rest, which ``at`` moves.
    """

    matrix: Array
    oxide: Array
    released: Array
    sidewall: Array
    gradient: Array
    blocks: tuple[HessianBlock, ...]
    dofs: tuple[Dof, ...]
    net: NDArray[np.intp]
    plates: tuple[Plate, ...]
    clusters: tuple[Cluster, ...]
    unassigned_area: float
    displacement: Array
    rest: RestGeometry

    def at(self, displacement: ArrayLike, *, progress: Progress | None = None) -> CapacitanceModel:
        """The model with its bodies displaced, equations (113) and (114).

        ``displacement`` holds one value for each degree of freedom in ``dofs``, measured
        from rest. Every matrix, the gradient, and the second derivative are those at the
        displaced bodies. Clusters, meshes, and windows stay as found at rest. Bodies that
        touch are an error. Stage sidewalls logs through the ``connflex_core`` logger and
        reports to ``progress`` as it does at rest.
        """
        moved = np.asarray(displacement, dtype=float)
        if moved.shape != (len(self.dofs),):
            raise ValueError(
                f"displacement must hold one value for each of {len(self.dofs)} degrees of freedom"
            )
        return _evaluate(self.rest, moved, progress)

    @property
    def net_matrix(self) -> Array:
        """Capacitance matrix between nets, with the electrodes of each net tied together."""
        count = int(np.max(self.net, initial=-1)) + 1
        tied = np.zeros((count, count))
        np.add.at(tied, (self.net[:, None], self.net[None, :]), self.matrix)
        return tied

    def force(self, potentials: ArrayLike) -> Array:
        """Generalized forces on ``dofs`` with ``potentials`` on the electrodes, equation (111)."""
        v = self._volts(potentials)
        force: Array = 0.5 * np.einsum("a,abk,b->k", v, self.gradient, v)
        return force

    def stiffness(self, potentials: ArrayLike) -> Array:
        """D by D electrostatic stiffness with ``potentials`` on the electrodes, equation (111).

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
        released, and sidewall parts and their total in farads, their nets, and the
        change of the total with the six degrees of freedom of each body that changes
        it. ``net_matrix``, ``clusters``, ``dofs``, ``displacement``, and
        ``unassigned_area`` are as on the model.
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
                }
                for cluster in self.clusters
            ],
            "dofs": [list(dof) for dof in self.dofs],
            "displacement": self.displacement.tolist(),
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
    """Capacitance between all electrodes of a structure at rest, equations (104) to (111).

    Plates between device and handle electrodes and sidewalls between device electrodes
    are each included unless ``params`` leaves them out. Sidewalls need
    ``sidewall_method`` and ``interaction_distance``, and ``element_size`` with method
    ``"field"``. Derivatives are in the six degrees of freedom of every movable body.
    Stages plates, clusters, and sidewalls log through the
    ``connflex_core`` logger and report to ``progress``; sidewalls advances once per
    cluster, and once per window for a cluster solved in windows.
    """
    params = params or ElectrostaticParams()
    plates = find_plates(structure, params, progress=progress) if params.include_plates else ()
    groups: list[RestCluster] = []
    if params.include_sidewalls:
        method, reach, size = (
            params.sidewall_method,
            params.interaction_distance,
            params.element_size,
        )
        if method is None or reach is None or (size is None and method == "field"):
            missing = (
                "sidewall_method"
                if method is None
                else "interaction_distance" if reach is None else "element_size"
            )
            raise ValueError(
                f"{missing} is required for sidewall capacitance; set it or include_sidewalls=False"
            )
        device = [e for e in structure.electrodes if e.layer == "device"]
        with stage("clusters", progress) as current:
            margin = structure.params.beam_overlap_margin or structure.params.heal_distance
            silicon = [_merged(electrode.region, margin) for electrode in device]
            graph: nx.Graph[int] = nx.Graph()
            graph.add_nodes_from(range(len(device)))
            near = STRtree(silicon).query(silicon, predicate="dwithin", distance=reach)
            graph.add_edges_from(near.T.tolist())
            for members in sorted(
                sorted(group) for group in nx.connected_components(graph) if len(group) > 1
            ):
                regions = tuple(silicon[k] for k in members)
                bodies = tuple(sorted({device[k].body for k in members}))
                mesh, cores = None, []
                if size is not None and method == "field":
                    mesh = mesh_outlines(regions, size, params.corner_grading)
                    if len(mesh.conductor) > params.max_elements:
                        cores = cut_cores(mesh, params.window_elements)
                groups.append(
                    RestCluster(
                        electrodes=tuple(device[k].index for k in members),
                        regions=regions,
                        bodies=bodies,
                        body=np.array([bodies.index(device[k].body) for k in members], np.intp),
                        mesh=mesh,
                        cores=cores,
                        windows=None,
                    )
                )
            current.count(clusters=len(groups), windows=sum(len(g.cores) for g in groups))
    rest = RestGeometry(structure, params, plates, tuple(groups))
    return _evaluate(
        rest, np.zeros(6 * sum(not body.anchor for body in structure.bodies)), progress
    )


def _evaluate(
    rest: RestGeometry, displacement: Array, progress: Progress | None
) -> CapacitanceModel:
    """Capacitance model of ``rest`` with its bodies displaced."""
    structure, params = rest.structure, rest.params
    count = len(structure.electrodes)
    oxide, released, sidewall = (np.zeros((count, count)) for _ in range(3))
    movable = [body.index for body in structure.bodies if not body.anchor]
    first = {index: 6 * position for position, index in enumerate(movable)}
    gradient = np.zeros((count, count, 6 * len(movable)))
    blocks: list[HessianBlock] = []

    def moved(body: int) -> Array:
        return displacement[first[body] : first[body] + 6] if body in first else np.zeros(6)

    for plate in rest.plates:
        value, slope, bend = plate.evaluate(moved(plate.bodies[0]), moved(plate.bodies[1]))
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
    steps = sum(max(len(group.cores), 1) for group in rest.groups)
    if rest.groups:
        with stage("sidewalls", progress, total=steps) as current:
            for group in rest.groups:
                centers = [structure.bodies[b].inertial.center for b in group.bodies]
                levers = [_in_plane(-thickness / 2 - center[2]) for center in centers]
                motion = np.array(
                    [lever @ moved(b) for lever, b in zip(levers, group.bodies, strict=True)]
                )
                pivots = np.array([center[:2] for center in centers]) + motion[:, :2]
                regions = _displaced(group, np.array([c[:2] for c in centers]), motion)
                if group.mesh is not None:
                    mesh = _carried(group, np.array([c[:2] for c in centers]), motion)
                    field, kept = solve_field(
                        mesh,
                        group.body,
                        pivots,
                        params,
                        quiet,
                        group.cores,
                        current.advance,
                        group.windows,
                    )
                    per_meter = (field.capacitance, field.gradient, field.hessian)
                    walls = FacingWalls.empty()
                    solved = Cluster(
                        group.electrodes, len(mesh.conductor), field.windows, field.asymmetry, walls
                    )
                    groups.append(replace(group, windows=kept) if group.cores else group)
                else:
                    between, slope, bend, walls = facing_plates(
                        regions, group.body, pivots, params.interaction_distance or 0.0
                    )
                    per_meter = (permittivity * between, permittivity * slope, permittivity * bend)
                    solved = Cluster(group.electrodes, len(walls.conductor), 0, 0.0, walls)
                    groups.append(group)
                clusters.append(solved)
                held, size = group.electrodes, len(group.electrodes)
                sidewall[np.ix_(held, held)] += thickness * per_meter[0]
                free = [k for k, index in enumerate(group.bodies) if index in first]
                if free:
                    travel = scipy.linalg.block_diag(*(levers[k] for k in free))
                    dofs = np.concatenate([first[group.bodies[k]] + np.arange(6) for k in free])
                    slope = per_meter[1][:, :, free].reshape(size, size, -1)
                    gradient[np.ix_(held, held, dofs)] += thickness * slope @ travel
                    bend = per_meter[2][:, :, free][:, :, :, :, free]
                    bend = bend.reshape(size, size, 3 * len(free), 3 * len(free))
                    values = thickness * np.einsum("abpq,pk,ql->abkl", bend, travel, travel)
                    blocks.append(HessianBlock(held, dofs, values))
                if not group.cores:
                    current.advance()

    bridges = {link.source for link in structure.links if link.kind == "bridge"}
    return CapacitanceModel(
        matrix=oxide + released + sidewall,
        oxide=oxide,
        released=released,
        sidewall=sidewall,
        gradient=gradient,
        blocks=tuple(blocks),
        dofs=tuple((index, name) for index in movable for name in DOF_NAMES),
        net=np.array([electrode.net for electrode in structure.electrodes], dtype=np.intp),
        plates=rest.plates,
        clusters=tuple(clusters),
        unassigned_area=sum(beam.polygon.area for beam in structure.beams)
        + sum(structure.rejected[index].polygon.area for index in bridges),
        displacement=displacement,
        rest=replace(rest, groups=tuple(groups)),
    )


def _turn(angle: float) -> Array:
    """Matrix that turns row vectors by ``angle``."""
    return np.array([[np.cos(angle), np.sin(angle)], [-np.sin(angle), np.cos(angle)]])


def _displaced(group: RestCluster, centers: Array, motion: Array) -> list[MultiPolygon]:
    """Silicon of a cluster with each body moved in the plane, equation (113)."""
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


def _carried(group: RestCluster, centers: Array, motion: Array) -> Mesh:
    """Mesh of a cluster with the elements of each body moved in the plane."""
    mesh = group.mesh
    if mesh is None:
        raise ValueError("cluster has no mesh")
    if not motion.any():
        return mesh
    of = group.body[mesh.conductor]
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


def _merged(region: MultiPolygon, margin: float) -> MultiPolygon:
    """Silicon of an electrode as one region, closed by ``margin`` where its polygons meet."""
    joined = shapely.union_all(shapely.get_parts(region))
    closed = joined.buffer(margin, join_style="mitre").buffer(-margin, join_style="mitre")
    return as_multipolygon(closed)


def _in_plane(lever: float) -> Array:
    """Travel along x and y and rotation about z at device mid-plane, equation (109)."""
    return np.array(
        [
            [1.0, 0.0, 0.0, 0.0, lever, 0.0],
            [0.0, 1.0, 0.0, -lever, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
        ]
    )


_BETWEEN = np.array([[1.0, -1.0], [-1.0, 1.0]])

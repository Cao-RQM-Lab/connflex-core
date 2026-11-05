"""Capacitance of a structure at rest."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import networkx as nx
import numpy as np
import scipy.linalg
import shapely
from numpy.typing import ArrayLike, NDArray
from shapely import MultiPolygon, STRtree

from connflex_core.electrostatics import VACUUM_PERMITTIVITY
from connflex_core.electrostatics.outlines import mesh_outlines
from connflex_core.electrostatics.plates import Plate, find_plates
from connflex_core.electrostatics.sidewalls import solve_mesh
from connflex_core.geometry.heal import as_multipolygon
from connflex_core.mechanics.assembly import DOF_NAMES, Dof
from connflex_core.params import ElectrostaticParams
from connflex_core.progress import Progress, quiet, stage

if TYPE_CHECKING:
    from connflex_core.model import Structure

Array = NDArray[np.float64]


@dataclass(frozen=True)
class Cluster:
    """Device electrodes whose sidewalls are solved as one field problem.

    Attributes:
        electrodes: Indices of its electrodes, ascending.
        elements: Number of boundary elements on their outlines.
        solved: False when the cluster was left out because it needs more than
            ``max_elements``.
    """

    electrodes: tuple[int, ...]
    elements: int
    solved: bool


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
class CapacitanceModel:
    """Capacitance between the electrodes of a structure at rest, in farads.

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
            and per solved cluster. Blocks add where they share entries.
        dofs: (body index, degree of freedom) for each of the D degrees of freedom: six
            per movable body in body order, as for the assembled mechanical matrices.
        net: Net of each electrode.
        plates: Plates behind ``oxide`` and ``released``.
        clusters: Groups of two or more device electrodes behind ``sidewall``.
        unassigned_area: Area in square meters of device silicon in beams and bridges,
            which belongs to no electrode and is not in the model.
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

    @property
    def net_matrix(self) -> Array:
        """Capacitance matrix between nets, with the electrodes of each net tied together."""
        count = int(np.max(self.net, initial=-1)) + 1
        tied = np.zeros((count, count))
        np.add.at(tied, (self.net[:, None], self.net[None, :]), self.matrix)
        return tied

    def force(self, potentials: ArrayLike) -> Array:
        """Generalized forces on ``dofs`` with ``potentials`` on the electrodes, equation (98)."""
        v = self._volts(potentials)
        force: Array = 0.5 * np.einsum("a,abk,b->k", v, self.gradient, v)
        return force

    def stiffness(self, potentials: ArrayLike) -> Array:
        """D by D electrostatic stiffness with ``potentials`` on the electrodes, equation (98).

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
        it. ``net_matrix``, ``clusters``, ``dofs``, and ``unassigned_area`` are as on
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
                    "solved": cluster.solved,
                }
                for cluster in self.clusters
            ],
            "dofs": [list(dof) for dof in self.dofs],
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
    """Capacitance between all electrodes of a structure, equations (91) to (98).

    Plates between device and handle electrodes and sidewalls between device electrodes
    are each included unless ``params`` leaves them out. Sidewalls need ``element_size``
    and ``interaction_distance``. Derivatives are in the six degrees of freedom of every
    movable body. Stages plates, clusters, and sidewalls log through the
    ``connflex_core`` logger and report to ``progress``; sidewalls advances once per
    cluster.
    """
    params = params or ElectrostaticParams()
    count = len(structure.electrodes)
    oxide, released, sidewall = (np.zeros((count, count)) for _ in range(3))
    movable = [body.index for body in structure.bodies if not body.anchor]
    first = {index: 6 * position for position, index in enumerate(movable)}
    gradient = np.zeros((count, count, 6 * len(movable)))
    blocks: list[HessianBlock] = []

    plates: tuple[Plate, ...] = ()
    if params.include_plates:
        plates = find_plates(structure, params, progress=progress)
        for plate in plates:
            pair = np.ix_(plate.electrodes, plate.electrodes)
            oxide[pair] += plate.oxide * _BETWEEN
            released[pair] += plate.released * _BETWEEN
            sides = [side for side, index in enumerate(plate.bodies) if index in first]
            if plate.bodies[0] == plate.bodies[1] or not sides:
                continue
            local = np.concatenate([np.arange(6 * side, 6 * side + 6) for side in sides])
            dofs = np.concatenate([first[plate.bodies[side]] + np.arange(6) for side in sides])
            _, slope, bend = plate.evaluate()
            gradient[np.ix_(plate.electrodes, plate.electrodes, dofs)] += (
                _BETWEEN[:, :, None] * slope[local]
            )
            values = _BETWEEN[:, :, None, None] * bend[np.ix_(local, local)]
            blocks.append(HessianBlock(plate.electrodes, dofs, values))

    clusters: list[Cluster] = []
    if params.include_sidewalls:
        size, reach = params.element_size, params.interaction_distance
        if size is None or reach is None:
            missing = "element_size" if size is None else "interaction_distance"
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
            groups = sorted(
                sorted(group) for group in nx.connected_components(graph) if len(group) > 1
            )
            current.count(clusters=len(groups))
        thickness = structure.stack.device.thickness
        midplane = -thickness / 2
        permittivity = VACUUM_PERMITTIVITY * params.ambient_permittivity
        with stage("sidewalls", progress, total=len(groups)) as current:
            for group in groups:
                electrodes = tuple(device[k].index for k in group)
                mesh = mesh_outlines([silicon[k] for k in group], size, params.corner_grading)
                solved = len(mesh.conductor) <= params.max_elements
                clusters.append(Cluster(electrodes, len(mesh.conductor), solved))
                if solved:
                    moving = sorted({device[k].body for k in group})
                    body = np.array([moving.index(device[k].body) for k in group], dtype=np.intp)
                    pivots = np.array([structure.bodies[b].inertial.center[:2] for b in moving])
                    field = solve_mesh(mesh, body, pivots, permittivity, quiet)
                    sidewall[np.ix_(electrodes, electrodes)] += thickness * field.capacitance
                    free = [k for k, index in enumerate(moving) if index in first]
                    if free:
                        lever = [
                            midplane - structure.bodies[moving[k]].inertial.center[2] for k in free
                        ]
                        travel = scipy.linalg.block_diag(*(_in_plane(arm) for arm in lever))
                        dofs = np.concatenate([first[moving[k]] + np.arange(6) for k in free])
                        slope = field.gradient[:, :, free].reshape(len(group), len(group), -1)
                        gradient[np.ix_(electrodes, electrodes, dofs)] += thickness * slope @ travel
                        bend = field.hessian[:, :, free][:, :, :, :, free]
                        bend = bend.reshape(len(group), len(group), 3 * len(free), 3 * len(free))
                        values = thickness * np.einsum("abpq,pk,ql->abkl", bend, travel, travel)
                        blocks.append(HessianBlock(electrodes, dofs, values))
                elif not params.skip_large_clusters:
                    raise ValueError(
                        f"electrodes {list(electrodes)} need {len(mesh.conductor)} boundary elements "
                        f"together, more than max_elements = {params.max_elements}: raise "
                        f"element_size or max_elements, lower interaction_distance, or set "
                        f"skip_large_clusters"
                    )
                current.advance()
            current.count(
                solved=sum(c.solved for c in clusters), skipped=sum(not c.solved for c in clusters)
            )

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
        plates=plates,
        clusters=tuple(clusters),
        unassigned_area=sum(beam.polygon.area for beam in structure.beams)
        + sum(structure.rejected[index].polygon.area for index in bridges),
    )


def _merged(region: MultiPolygon, margin: float) -> MultiPolygon:
    """Silicon of an electrode as one region, closed by ``margin`` where its polygons meet."""
    joined = shapely.union_all(shapely.get_parts(region))
    closed = joined.buffer(margin, join_style="mitre").buffer(-margin, join_style="mitre")
    return as_multipolygon(closed)


def _in_plane(lever: float) -> Array:
    """Travel along x and y and rotation about z at device mid-plane, equation (96)."""
    return np.array(
        [
            [1.0, 0.0, 0.0, 0.0, lever, 0.0],
            [0.0, 1.0, 0.0, -lever, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
        ]
    )


_BETWEEN = np.array([[1.0, -1.0], [-1.0, 1.0]])

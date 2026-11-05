"""Capacitance of a structure at rest."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import networkx as nx
import numpy as np
import shapely
from numpy.typing import NDArray
from shapely import MultiPolygon, STRtree

from connflex_core.electrostatics import VACUUM_PERMITTIVITY
from connflex_core.electrostatics.outlines import mesh_outlines
from connflex_core.electrostatics.plates import Plate, find_plates
from connflex_core.electrostatics.sidewalls import solve_mesh
from connflex_core.geometry.heal import as_multipolygon
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
        net: Net of each electrode.
        plates: Plates behind ``oxide`` and ``released``.
        clusters: Groups of two or more device electrodes behind ``sidewall``.
    """

    matrix: Array
    oxide: Array
    released: Array
    sidewall: Array
    net: NDArray[np.intp]
    plates: tuple[Plate, ...]
    clusters: tuple[Cluster, ...]

    @property
    def net_matrix(self) -> Array:
        """Capacitance matrix between nets, with the electrodes of each net tied together."""
        count = int(np.max(self.net, initial=-1)) + 1
        tied = np.zeros((count, count))
        np.add.at(tied, (self.net[:, None], self.net[None, :]), self.matrix)
        return tied


def find_capacitance(
    structure: Structure,
    params: ElectrostaticParams | None = None,
    *,
    progress: Progress | None = None,
) -> CapacitanceModel:
    """Capacitance between all electrodes of a structure, equations (91) to (95).

    Plates between device and handle electrodes and sidewalls between device electrodes
    are each included unless ``params`` leaves them out. Sidewalls need ``element_size``
    and ``interaction_distance``. Stages plates, clusters, and sidewalls log through the
    ``connflex_core`` logger and report to ``progress``; sidewalls advances once per
    cluster.
    """
    params = params or ElectrostaticParams()
    count = len(structure.electrodes)
    oxide, released, sidewall = (np.zeros((count, count)) for _ in range(3))

    plates: tuple[Plate, ...] = ()
    if params.include_plates:
        plates = find_plates(structure, params, progress=progress)
        for plate in plates:
            pair = np.ix_(plate.electrodes, plate.electrodes)
            oxide[pair] += plate.oxide * _BETWEEN
            released[pair] += plate.released * _BETWEEN

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

    return CapacitanceModel(
        matrix=oxide + released + sidewall,
        oxide=oxide,
        released=released,
        sidewall=sidewall,
        net=np.array([electrode.net for electrode in structure.electrodes], dtype=np.intp),
        plates=plates,
        clusters=tuple(clusters),
    )


def _merged(region: MultiPolygon, margin: float) -> MultiPolygon:
    """Silicon of an electrode as one region, closed by ``margin`` where its polygons meet."""
    joined = shapely.union_all(shapely.get_parts(region))
    closed = joined.buffer(margin, join_style="mitre").buffer(-margin, join_style="mitre")
    return as_multipolygon(closed)


_BETWEEN = np.array([[1.0, -1.0], [-1.0, 1.0]])

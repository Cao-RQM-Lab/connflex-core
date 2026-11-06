"""Nodes of the electrical network that sources and signals see."""

from __future__ import annotations

from collections.abc import Collection, Sequence
from itertools import pairwise
from typing import TYPE_CHECKING

import networkx as nx
import numpy as np
from numpy.typing import NDArray

from connflex_core.electrical.electrodes import Electrode

if TYPE_CHECKING:
    from connflex_core.model import Structure

Array = NDArray[np.float64]


def electrode_indices(
    structure: Structure, given: Sequence[object], floating: Collection[int], role: str
) -> list[int]:
    """Indices of electrodes given as objects of ``structure`` on nets that do not float."""
    count = len(structure.electrodes)
    places = []
    for item in given:
        if not isinstance(item, Electrode):
            raise TypeError(
                f"a {role} is an Electrode, such as structure.electrode_at(point, layer), not "
                f"{type(item).__name__} {item!r}: electrodes and nets are numbered apart"
            )
        if not (0 <= item.index < count and structure.electrodes[item.index] is item):
            raise ValueError(f"electrode {item.index} is of another structure")
        if item.net in floating:
            raise ValueError(
                f"electrode {item.index} is on floating net {item.net}, which no source holds"
            )
        places.append(item.index)
    return places


def network_nodes(
    structure: Structure, floating: Collection[int], ports: Sequence[int], role: str
) -> tuple[list[tuple[int, ...]], list[tuple[int, ...]], Array]:
    """Nodes of nets with ports, ports first, nodes of floating nets, and conductance
    between the first, equation (130).

    Links without resistance join electrodes into one node, and every floating net is
    one node. Electrodes of other nets belong to no node.
    """
    net = np.array([electrode.net for electrode in structure.electrodes], dtype=np.intp)
    sourced = {int(net[port]) for port in ports}
    shorts: nx.Graph[int] = nx.Graph()
    shorts.add_nodes_from(range(len(net)))
    for link in structure.links:
        if not link.resistance:
            shorts.add_edge(*link.electrodes)
    for each in floating:
        shorts.add_edges_from(pairwise(np.flatnonzero(net == each).tolist()))
    groups = {a: tuple(sorted(group)) for group in nx.connected_components(shorts) for a in group}
    nodes: list[tuple[int, ...]] = []
    for port in ports:
        if groups[port] in nodes:
            raise ValueError(
                f"{role}s at electrodes {nodes[nodes.index(groups[port])][0]} and {port} are "
                f"one node: links without resistance join them"
            )
        nodes.append(groups[port])
    inner = {groups[a] for a in range(len(net)) if net[a] in sourced and groups[a] not in nodes}
    nodes.extend(sorted(inner))
    hidden = sorted({groups[a] for a in range(len(net)) if net[a] in floating})
    index = {a: k for k, node in enumerate(nodes) for a in node}
    conductance = np.zeros((len(nodes), len(nodes)))
    for link in structure.links:
        a, b = link.electrodes
        if link.resistance and a in index and index[a] != index[b]:
            pair = [index[a], index[b]]
            conductance[np.ix_(pair, pair)] += np.array([[1, -1], [-1, 1]]) / link.resistance
    return nodes, hidden, conductance

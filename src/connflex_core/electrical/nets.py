"""Nets from conductors and vias."""

from __future__ import annotations

import warnings
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any, Literal

import networkx as nx
from shapely import Point, Polygon, STRtree
from shapely.geometry import mapping

from connflex_core.geometry.heal import overlapping


@dataclass(frozen=True)
class Conductor:
    """A connected region of one conducting layer.

    Attributes:
        index: Position in ``Nets.conductors``.
        layer: ``"device"`` or ``"handle"``.
        polygon: The region, in meters.
        net: Index of the conductor's net.
    """

    index: int
    layer: str
    polygon: Polygon
    net: int


@dataclass(frozen=True)
class Via:
    """A via and the conductors it overlaps on each side.

    Attributes:
        index: Position in ``Nets.vias``.
        polygon: The via region, in meters.
        device: Indices of device conductors it overlaps.
        handle: Indices of handle conductors it overlaps. A via with no conductor on
            one side joins nothing.
    """

    index: int
    polygon: Polygon
    device: tuple[int, ...]
    handle: tuple[int, ...]


LAYERS = ("device", "handle")
End = tuple[float, float] | tuple[float, float, Literal["device", "handle"]]


@dataclass(frozen=True)
class Wire:
    """Connection made off the die between two pads, such as a wire bond or a trace.

    Attributes:
        start: One end: a point (x, y) in meters on a pad of the device layer, or
            (x, y, layer) with layer ``"device"`` or ``"handle"``.
        end: The other end, in the same form.
        resistance: Resistance between the two ends, in ohms. Zero is ideal.
    """

    start: End
    end: End
    resistance: float = 0.0

    def __post_init__(self) -> None:
        for end in (self.start, self.end):
            if len(end) not in (2, 3) or (len(end) == 3 and end[2] not in LAYERS):
                raise ValueError(
                    f"wire end must be (x, y) or (x, y, layer) with layer in {LAYERS}, got {end}"
                )
        if not self.resistance >= 0:
            raise ValueError(f"resistance must be at least 0, got {self.resistance}")

    @property
    def ends(self) -> tuple[tuple[Point, str], tuple[Point, str]]:
        """Point and layer of each end."""
        start, end = (
            (Point(e[0], e[1]), "device" if len(e) == 2 else str(e[2]))
            for e in (self.start, self.end)
        )
        return start, end


@dataclass(frozen=True)
class Connection:
    """A wire and the conductors under its two ends.

    Attributes:
        index: Position in ``Nets.wires``, the same as in the wires passed to extraction.
        wire: The wire.
        conductors: Indices of the conductors under its start and its end.
    """

    index: int
    wire: Wire
    conductors: tuple[int, int]


@dataclass(frozen=True)
class Nets:
    """Electrical nets: conductors joined by vias and by wires off the die.

    Attributes:
        conductors: Device conductors first, then handle conductors.
        vias: All vias, in input order.
        count: Number of nets.
        wires: Every wire, in input order, with the conductors at its ends.
    """

    conductors: tuple[Conductor, ...]
    vias: tuple[Via, ...]
    count: int
    wires: tuple[Connection, ...] = ()

    @property
    def graph(self) -> nx.MultiGraph[int]:
        """Conductors as nodes; vias as edges between each device and handle conductor they
        join, keyed by via index, and wires as edges between the conductors at their
        ends, keyed by ``("wire", index)``."""
        graph: nx.MultiGraph[int] = nx.MultiGraph()
        for conductor in self.conductors:
            graph.add_node(conductor.index, conductor=conductor)
        for via in self.vias:
            for a in via.device:
                for b in via.handle:
                    graph.add_edge(a, b, key=via.index, via=via)
        for connection in self.wires:
            graph.add_edge(*connection.conductors, key=("wire", connection.index), wire=connection)
        return graph

    def conductor_at(self, point: Point | tuple[float, float], layer: str) -> Conductor:
        """Conductor of ``layer`` covering ``point``."""
        point = Point(point)
        for conductor in self.conductors:
            if conductor.layer == layer and conductor.polygon.covers(point):
                return conductor
        raise ValueError(f"no {layer} conductor at ({point.x}, {point.y})")

    def to_dict(self) -> dict[str, Any]:
        """JSON-serializable mapping with polygons as GeoJSON-style mappings."""
        return {
            "count": self.count,
            "conductors": [
                {
                    "index": conductor.index,
                    "layer": conductor.layer,
                    "net": conductor.net,
                    "polygon": mapping(conductor.polygon),
                }
                for conductor in self.conductors
            ],
            "vias": [
                {
                    "index": via.index,
                    "polygon": mapping(via.polygon),
                    "device": list(via.device),
                    "handle": list(via.handle),
                }
                for via in self.vias
            ],
            "wires": [
                {
                    "index": connection.index,
                    "start": list(connection.wire.start),
                    "end": list(connection.wire.end),
                    "resistance": connection.wire.resistance,
                    "conductors": list(connection.conductors),
                }
                for connection in self.wires
            ],
        }


def find_nets(
    device: Sequence[Polygon],
    handle: Sequence[Polygon],
    vias: Sequence[Polygon] = (),
    wires: Sequence[Wire] = (),
) -> Nets:
    """Group device and handle conductors into nets through vias and wires.

    Each polygon is one conductor. A via joins every device conductor and every handle
    conductor it overlaps; a via overlapping no conductor on one side is ignored with a warning.
    A wire joins the conductors that cover its two ends; an end that no conductor of its
    layer covers is an error. Nets are numbered in order of their first conductor.
    """
    parent = list(range(len(device) + len(handle)))

    def root(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    trees = (STRtree(list(device)), STRtree(list(handle)))
    joined: list[Via] = []
    for index, via in enumerate(vias):
        on_device = _overlapping(trees[0], device, via, 0)
        on_handle = _overlapping(trees[1], handle, via, len(device))
        joined.append(Via(index, via, on_device, on_handle))
        if not on_device or not on_handle:
            missing = "device" if not on_device else "handle"
            warnings.warn(f"via {index} overlaps no {missing} conductor; ignored", stacklevel=2)
            continue
        members = [*on_device, *on_handle]
        first = root(members[0])
        for member in members[1:]:
            parent[root(member)] = first

    connections: list[Connection] = []
    for index, wire in enumerate(wires):
        under = []
        for point, layer in wire.ends:
            tree, offset = (trees[0], 0) if layer == "device" else (trees[1], len(device))
            polygons = device if layer == "device" else handle
            hits = [int(j) for j in tree.query(point) if polygons[j].covers(point)]
            if not hits:
                raise ValueError(
                    f"wire {index}: end at ({point.x:g}, {point.y:g}) m is on no {layer} "
                    f"conductor: empty space, or silicon dropped as unconnected to the anchor"
                )
            under.append(offset + hits[0])
        connections.append(Connection(index, wire, (under[0], under[1])))
        parent[root(under[1])] = root(under[0])

    numbers: dict[int, int] = {}
    conductors: list[Conductor] = []
    layered = [*(("device", p) for p in device), *(("handle", p) for p in handle)]
    for index, (layer, polygon) in enumerate(layered):
        net = numbers.setdefault(root(index), len(numbers))
        conductors.append(Conductor(index, layer, polygon, net))
    return Nets(tuple(conductors), tuple(joined), len(numbers), tuple(connections))


def _overlapping(
    tree: STRtree, polygons: Sequence[Polygon], via: Polygon, offset: int
) -> tuple[int, ...]:
    return tuple(sorted(offset + j for j in overlapping(tree, polygons, via)))

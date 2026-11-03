"""Nets from conductors and vias."""

from __future__ import annotations

import warnings
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import networkx as nx
from shapely import Point, Polygon, STRtree
from shapely.geometry import mapping


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


@dataclass(frozen=True)
class Nets:
    """Electrical nets: conductors joined by vias.

    Attributes:
        conductors: Device conductors first, then handle conductors.
        vias: All vias, in input order.
        count: Number of nets.
    """

    conductors: tuple[Conductor, ...]
    vias: tuple[Via, ...]
    count: int

    @property
    def graph(self) -> nx.MultiGraph[int]:
        """Conductors as nodes and vias as edges between each device and handle conductor they join."""
        graph: nx.MultiGraph[int] = nx.MultiGraph()
        for conductor in self.conductors:
            graph.add_node(conductor.index, conductor=conductor)
        for via in self.vias:
            for a in via.device:
                for b in via.handle:
                    graph.add_edge(a, b, key=via.index, via=via)
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
        }


def find_nets(
    device: Sequence[Polygon], handle: Sequence[Polygon], vias: Sequence[Polygon] = ()
) -> Nets:
    """Group device and handle conductors into nets through the vias that overlap them.

    Each polygon is one conductor. A via joins every device conductor and every handle
    conductor it overlaps; a via overlapping no conductor on one side is ignored with a warning.
    Nets are numbered in order of their first conductor.
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

    numbers: dict[int, int] = {}
    conductors: list[Conductor] = []
    layered = [*(("device", p) for p in device), *(("handle", p) for p in handle)]
    for index, (layer, polygon) in enumerate(layered):
        net = numbers.setdefault(root(index), len(numbers))
        conductors.append(Conductor(index, layer, polygon, net))
    return Nets(tuple(conductors), tuple(joined), len(numbers))


def _overlapping(
    tree: STRtree, polygons: Sequence[Polygon], via: Polygon, offset: int
) -> tuple[int, ...]:
    hits = tree.query(via, predicate="intersects")
    return tuple(sorted(offset + int(j) for j in hits if via.intersection(polygons[j]).area > 0))

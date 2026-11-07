"""Anchor selection and connectivity over the layer stack."""

from __future__ import annotations

from collections import deque
from collections.abc import Callable, Iterator, Sequence

from shapely import Point, Polygon, STRtree

from connflex_core.geometry.heal import overlapping

Node = tuple[int, int]


def select_anchor(handle: Sequence[Polygon], anchor: int | tuple[float, float] | None) -> int:
    """Index of the anchor handle polygon: the largest, the given index, or the one containing a point."""
    if anchor is None:
        return max(range(len(handle)), key=lambda index: handle[index].area)
    if isinstance(anchor, int):
        if not 0 <= anchor < len(handle):
            raise ValueError(
                f"anchor index {anchor} is out of range; the handle layer has {len(handle)} polygons"
            )
        return anchor
    point = Point(anchor)
    for index, polygon in enumerate(handle):
        if polygon.contains(point):
            return index
    raise ValueError(f"anchor point {anchor} lies in no handle polygon")


def layered_components(
    layers: Sequence[Sequence[Polygon]], seed: Node, advance: Callable[[], None] | None = None
) -> list[list[Node]]:
    """Connected components where polygons in adjacent layers connect when they overlap.

    Nodes are ``(layer, index)`` pairs. The component containing ``seed`` comes first;
    others follow in layer and index order. ``advance`` is called once per polygon claimed.
    """
    trees = [STRtree(list(layer)) for layer in layers]
    visited: set[Node] = set()
    components: list[list[Node]] = []
    starts = [seed, *((layer, i) for layer, parts in enumerate(layers) for i in range(len(parts)))]
    for start in starts:
        if start in visited:
            continue
        component: list[Node] = []
        queue = deque([start])
        visited.add(start)
        while queue:
            node = queue.popleft()
            component.append(node)
            if advance is not None:
                advance()
            for neighbor in _overlapping_neighbors(layers, trees, node):
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(neighbor)
        components.append(component)
    return components


def _overlapping_neighbors(
    layers: Sequence[Sequence[Polygon]], trees: Sequence[STRtree], node: Node
) -> Iterator[Node]:
    layer, index = node
    polygon = layers[layer][index]
    for adjacent in (layer - 1, layer + 1):
        if 0 <= adjacent < len(layers):
            for j in overlapping(trees[adjacent], layers[adjacent], polygon):
                yield (adjacent, j)

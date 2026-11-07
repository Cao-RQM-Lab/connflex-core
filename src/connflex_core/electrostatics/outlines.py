"""Boundary elements on conductor outlines."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import shapely
from numpy.typing import NDArray
from shapely import MultiPolygon, Polygon, STRtree

Array = NDArray[np.float64]


@dataclass(frozen=True, eq=False)
class Mesh:
    """Straight boundary elements on conductor outlines.

    Attributes:
        start: N by 2 start points of the elements, in meters.
        end: N by 2 end points.
        conductor: Index of the conductor each element lies on.
    """

    start: Array
    end: Array
    conductor: NDArray[np.intp]

    @classmethod
    def empty(cls) -> Mesh:
        return cls(np.zeros((0, 2)), np.zeros((0, 2)), np.zeros(0, dtype=np.intp))


VANISHING = 1e-6


def mesh_outlines(conductors: Sequence[MultiPolygon], size: float, grading: float) -> Mesh:
    """Divide outlines of ``conductors`` into elements, equations (72) to (74).

    Elements are at most ``size`` long and shorten toward convex corners, down to
    ``size / grading`` at a right angle or sharper. A vertex within ``VANISHING`` times
    ``size`` of the next, as rounding leaves them, is dropped. A hole is meshed only
    when it holds another polygon: an empty cavity in a conductor holds no field.
    """
    polygons = [polygon for conductor in conductors for polygon in conductor.geoms]
    owner = [index for index, conductor in enumerate(conductors) for _ in conductor.geoms]
    tree = STRtree(polygons)
    start: list[NDArray[np.complex128]] = []
    end: list[NDArray[np.complex128]] = []
    conductor: list[int] = []
    for position, polygon in enumerate(polygons):
        rings = [(polygon.exterior, 1.0)]
        for ring in polygon.interiors:
            inside = tree.query(Polygon(ring), predicate="intersects")
            if np.any(inside != position):
                rings.append((ring, -1.0))
        for ring, role in rings:
            for points in _ring_nodes(shapely.get_coordinates(ring), role, size, grading):
                start.append(points[:-1])
                end.append(points[1:])
                conductor.extend([owner[position]] * (points.size - 1))
    first, last = np.concatenate(start), np.concatenate(end)
    return Mesh(
        start=np.column_stack([first.real, first.imag]),
        end=np.column_stack([last.real, last.imag]),
        conductor=np.array(conductor, dtype=np.intp),
    )


def _ring_nodes(
    coordinates: Array, role: float, size: float, grading: float
) -> list[NDArray[np.complex128]]:
    """Nodes along each edge of a closed ring, as complex numbers x + iy.

    ``role`` is 1 for an exterior ring and -1 for a hole, where silicon lies outside.
    """
    vertex = coordinates[:-1, 0] + 1j * coordinates[:-1, 1]
    vertex = vertex[np.abs(np.roll(vertex, -1) - vertex) > VANISHING * size]
    edge = np.roll(vertex, -1) - vertex
    twice_area = np.sum(
        vertex.real * np.roll(vertex.imag, -1) - np.roll(vertex.real, -1) * vertex.imag
    )
    turn = np.angle(edge / np.roll(edge, 1))
    convex = np.clip(role * np.sign(twice_area) * turn / (np.pi / 2), 0.0, 1.0)
    corner = size / grading**convex
    return [
        vertex[k]
        + edge[k] * edge_fractions(abs(edge[k]), corner[k], corner[k - len(vertex) + 1], size)
        for k in range(len(vertex))
    ]


def edge_fractions(length: float, left: float, right: float, size: float) -> Array:
    """Node positions as fractions of an edge whose corners ask for ``left`` and ``right``.

    Target element length at distance x from the left corner is
    min(size, left + x, right + length - x); nodes sit at equal steps of the integral of
    its inverse, equations (73) and (74).
    """
    rise = min(max(size - left, 0.0), length)
    fall = max(min(length - (size - right), length), 0.0)
    if rise > fall:
        rise = fall = min(max((length + right - left) / 2, 0.0), length)
    climbing = np.log1p(rise / left)
    level = (fall - rise) / size
    falling = np.log1p((length - fall) / right)
    total = climbing + level + falling
    count = max(1, int(np.ceil(total - 1e-9)))
    step = total * np.arange(count + 1) / count
    position = np.where(
        step <= climbing,
        left * np.expm1(np.minimum(step, climbing)),
        np.where(
            step <= climbing + level,
            rise + (step - climbing) * size,
            length - right * np.expm1(np.clip(total - step, 0.0, falling)),
        ),
    )
    position[0], position[-1] = 0.0, length
    fractions: Array = position / length
    return fractions

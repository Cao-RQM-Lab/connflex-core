"""Morphological operations on polygon layers."""

from __future__ import annotations

import shapely
from shapely import MultiPolygon, Polygon
from shapely.geometry.base import BaseGeometry


def as_multipolygon(geometry: BaseGeometry) -> MultiPolygon:
    """Collect the polygonal parts of any geometry into a MultiPolygon."""
    parts = [
        part
        for part in shapely.get_parts(geometry)
        if isinstance(part, Polygon) and not part.is_empty
    ]
    return MultiPolygon(parts)


def heal(geometry: MultiPolygon, distance: float, tolerance: float) -> MultiPolygon:
    """Close gaps narrower than twice ``distance``, then simplify with ``tolerance``."""
    closed = geometry.buffer(distance, join_style="mitre").buffer(-distance, join_style="mitre")
    return as_multipolygon(closed.simplify(tolerance, preserve_topology=True))


def opening(geometry: MultiPolygon, radius: float) -> MultiPolygon:
    """Remove features narrower than twice ``radius``; wider features come back unchanged."""
    eroded = geometry.buffer(-radius, join_style="mitre")
    return as_multipolygon(eroded.buffer(radius, join_style="mitre"))


def snap(geometry: MultiPolygon, grid: float) -> MultiPolygon:
    """Snap coordinates to a grid."""
    return as_multipolygon(shapely.set_precision(geometry, grid))

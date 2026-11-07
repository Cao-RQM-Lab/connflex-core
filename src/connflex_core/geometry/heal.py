"""Morphological operations on polygon layers."""

from __future__ import annotations

import numpy as np
import shapely
from numpy.typing import NDArray
from shapely import MultiPolygon, Polygon
from shapely.geometry.base import BaseGeometry

RESIDUE = 1e-12


def as_multipolygon(geometry: BaseGeometry | NDArray[np.object_]) -> MultiPolygon:
    """Collect polygonal parts of a geometry, or of an array of geometries, into a MultiPolygon."""
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
    """Remove features narrower than twice ``radius``; wider features come back unchanged.

    Each polygon is opened on its own. The opening of a polygon lies inside it, so the parts
    of a disjoint set stay disjoint and the union that one buffer over the whole set would
    perform is skipped. Vertices within ``RESIDUE`` of the line through their neighbors,
    which the mitre buffers leave along edges at an angle to the axes where a feature was
    removed, are dropped.
    """
    eroded = shapely.buffer(shapely.get_parts(geometry), -radius, join_style="mitre")
    dilated = shapely.buffer(eroded, radius, join_style="mitre")
    return as_multipolygon(shapely.simplify(dilated, RESIDUE, preserve_topology=True))


def snap(geometry: MultiPolygon, grid: float) -> MultiPolygon:
    """Snap coordinates to a grid."""
    return as_multipolygon(shapely.set_precision(geometry, grid))

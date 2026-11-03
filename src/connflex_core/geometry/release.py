"""Release holes and surviving buried oxide."""

from __future__ import annotations

from collections.abc import Sequence

import shapely
from shapely import MultiPolygon, Polygon

from connflex_core.geometry.heal import as_multipolygon


def fill_release_holes(device: MultiPolygon, max_area: float) -> tuple[MultiPolygon, list[Polygon]]:
    """Fill holes of at most ``max_area`` in device polygons and return them separately."""
    filled: list[Polygon] = []
    holes: list[Polygon] = []
    for polygon in device.geoms:
        kept = []
        for ring in polygon.interiors:
            hole = Polygon(ring)
            if hole.area <= max_area:
                holes.append(hole)
            else:
                kept.append(ring)
        filled.append(Polygon(polygon.exterior, kept))
    return MultiPolygon(filled), holes


def surviving_box(
    device: MultiPolygon,
    handle: MultiPolygon,
    release: MultiPolygon | None,
    holes: Sequence[Polygon],
    distance: float,
    resolution: int,
) -> MultiPolygon:
    """Oxide protected by both layers, minus openings, eroded by undercut ``distance``."""
    protected = device.intersection(handle)
    openings = [*holes, *(release.geoms if release is not None else [])]
    if openings:
        protected = protected.difference(shapely.union_all(openings))
    return as_multipolygon(protected.buffer(-distance, quad_segs=resolution))

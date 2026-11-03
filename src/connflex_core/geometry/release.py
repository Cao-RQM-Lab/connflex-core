"""Release holes and surviving buried oxide."""

from __future__ import annotations

from shapely import MultiPolygon, Polygon

from connflex_core.geometry.heal import as_multipolygon


def fill_release_holes(device: MultiPolygon, max_area: float) -> tuple[MultiPolygon, list[Polygon]]:
    """Fill holes of at most ``max_area`` in device polygons and return them separately.

    Called on the raw layout, so healing runs on the filled outlines rather than on every
    hole outline, and again on the healed layout for cavities that healing closes.
    """
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
    openings: MultiPolygon,
    distance: float,
    resolution: int,
) -> MultiPolygon:
    """Oxide protected by both layers, minus ``openings``, eroded by undercut ``distance``.

    ``openings`` holds the release holes and the release layer as disjoint polygons. It is
    built once by ``extract`` and shared with the mass properties.
    """
    protected = device.intersection(handle)
    if not openings.is_empty:
        protected = protected.difference(openings)
    return as_multipolygon(protected.buffer(-distance, quad_segs=resolution))

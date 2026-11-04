"""Release holes and surviving buried oxide."""

from __future__ import annotations

from collections.abc import Sequence

import shapely
from shapely import MultiPolygon, Polygon, STRtree

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

    Eroding a difference equals the eroded minuend minus the dilated subtrahend, so the
    protected region is eroded on its own and only the openings within ``distance`` of the
    result are dilated and subtracted, instead of eroding one polygon that carries every
    opening as a hole.
    """
    eroded = device.intersection(handle).buffer(-distance, quad_segs=resolution)
    if openings.is_empty or eroded.is_empty:
        return as_multipolygon(eroded)
    parts = shapely.get_parts(openings)
    near = STRtree(parts).query(eroded, predicate="dwithin", distance=distance)
    if near.size == 0:
        return as_multipolygon(eroded)
    dilated = shapely.buffer(parts[near], distance, quad_segs=resolution)
    return as_multipolygon(eroded.difference(shapely.union_all(dilated)))


def openings_within(
    tree: STRtree, openings: Sequence[Polygon], geometry: MultiPolygon
) -> list[Polygon]:
    """Openings inside ``geometry``: those it covers as drawn, those crossing its edge clipped.

    ``tree`` indexes ``openings``. The result is disjoint, so its area moments can be
    subtracted from those of ``geometry``.
    """
    inside: list[Polygon] = []
    for part in geometry.geoms:
        covered = {int(i) for i in tree.query(part, predicate="covers")}
        inside.extend(openings[i] for i in sorted(covered))
        for i in tree.query(part, predicate="intersects"):
            if int(i) not in covered:
                inside.extend(as_multipolygon(openings[int(i)].intersection(part)).geoms)
    return inside

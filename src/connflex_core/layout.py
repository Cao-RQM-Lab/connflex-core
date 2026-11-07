"""Canonical layout: polygon layers in meters."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

import numpy as np
import shapely
from numpy.typing import NDArray
from shapely import MultiPolygon, Polygon, STRtree
from shapely.validation import explain_validity

from connflex_core.geometry.heal import polygon_parts

ROLES = ("device", "handle", "release", "vias")

PolygonInput = Polygon | MultiPolygon | Iterable[Polygon | MultiPolygon]


@dataclass(frozen=True)
class LayerMap:
    """Layer and datatype of each layout role, for GDS, KLayout, and gdsfactory readers.

    Attributes:
        device: Device layer and datatype.
        handle: Handle layer and datatype.
        release: Release-opening layer and datatype, or None when there are none.
        vias: Via layer and datatype, or None when there are none.
    """

    device: tuple[int, int]
    handle: tuple[int, int]
    release: tuple[int, int] | None = None
    vias: tuple[int, int] | None = None


@dataclass(frozen=True)
class Layout:
    """Polygon layers of one layout, in meters.

    Every polygon marks where the feature named by its role is present. Polygons within one
    role that overlap or touch are merged; the rest keep their order. Each polygon is stored
    in canonical vertex order, so results depend on shapes alone and not on where a ring
    starts, which way it runs, or how its holes are listed.

    Attributes:
        device: Device-layer silicon.
        handle: Handle-layer silicon.
        release: Release openings, or None.
        vias: Regions where device and handle layers are electrically connected, or None.
    """

    device: MultiPolygon
    handle: MultiPolygon
    release: MultiPolygon | None = None
    vias: MultiPolygon | None = None

    @classmethod
    def from_shapely(
        cls,
        *,
        device: PolygonInput,
        handle: PolygonInput,
        release: PolygonInput | None = None,
        vias: PolygonInput | None = None,
        unit: float,
    ) -> Layout:
        """Build a layout from shapely polygons drawn in ``unit`` meters per coordinate unit."""
        from connflex_core.adapters.shapely_ import from_shapely

        roles: dict[str, PolygonInput] = {"device": device, "handle": handle}
        if release is not None:
            roles["release"] = release
        if vias is not None:
            roles["vias"] = vias
        return from_shapely(roles, unit=unit)


def normalize(role: str, polygons: PolygonInput, unit: float) -> MultiPolygon:
    """Validate, merge, scale, and canonicalize one role's polygons.

    Polygons that overlap or touch are unioned. The others are kept as given, in order, so a
    role of many separate polygons costs one spatial query instead of a union. Every polygon
    is then normalized, since simplification later in the pipeline depends on where a ring
    starts.
    """
    if unit <= 0:
        raise ValueError(f"unit must be positive, got {unit}")
    items = [polygons] if isinstance(polygons, Polygon | MultiPolygon) else list(polygons)
    for index in np.flatnonzero(~shapely.is_valid(np.array(items, dtype=object))):
        raise ValueError(f"{role} polygon {index} is invalid: {explain_validity(items[index])}")
    parts = np.array(polygon_parts(np.array(items, dtype=object)), dtype=object)
    if parts.size == 0:
        raise ValueError(f"{role} has no polygons")
    scaled = shapely.transform(_merge_touching(parts), lambda coords: coords * unit)
    return MultiPolygon(list(shapely.normalize(scaled)))


def _merge_touching(parts: NDArray[np.object_]) -> NDArray[np.object_]:
    left, right = STRtree(parts).query(parts, predicate="intersects")
    touching = np.unique(left[left != right])
    if touching.size == 0:
        return parts
    keep = np.ones(parts.size, dtype=bool)
    keep[touching] = False
    merged = np.array(polygon_parts(shapely.union_all(parts[touching])), dtype=object)
    return np.concatenate([parts[keep], merged])

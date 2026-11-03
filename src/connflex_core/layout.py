"""Canonical layout: polygon layers in meters."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

import shapely
from shapely import MultiPolygon, Polygon
from shapely.validation import explain_validity

ROLES = ("device", "handle", "release", "vias")

PolygonInput = Polygon | MultiPolygon | Iterable[Polygon | MultiPolygon]


@dataclass(frozen=True)
class LayerMap:
    """GDS layer and datatype of each layout role.

    Attributes:
        device: Layer and datatype of the device layer.
        handle: Layer and datatype of the handle layer.
        release: Layer and datatype of the release openings, or None when there are none.
        vias: Layer and datatype of the vias, or None when there are none.
    """

    device: tuple[int, int]
    handle: tuple[int, int]
    release: tuple[int, int] | None = None
    vias: tuple[int, int] | None = None


@dataclass(frozen=True)
class Layout:
    """Polygon layers of one layout, in meters.

    Every polygon marks where the feature named by its role is present. Polygons within one
    role are unioned.

    Attributes:
        device: Device-layer silicon.
        handle: Handle-layer silicon.
        release: Release openings, or None.
        vias: Regions where the device layer is electrically connected to the handle, or None.
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
    """Validate, union, and scale the polygons of one role."""
    if unit <= 0:
        raise ValueError(f"unit must be positive, got {unit}")
    items = [polygons] if isinstance(polygons, Polygon | MultiPolygon) else list(polygons)
    for index, geometry in enumerate(items):
        if not geometry.is_valid:
            raise ValueError(f"{role} polygon {index} is invalid: {explain_validity(geometry)}")
    union = shapely.transform(shapely.union_all(items), lambda coords: coords * unit)
    parts = [part for part in shapely.get_parts(union) if isinstance(part, Polygon)]
    if not parts or all(part.is_empty for part in parts):
        raise ValueError(f"{role} has no polygons")
    return MultiPolygon(parts)

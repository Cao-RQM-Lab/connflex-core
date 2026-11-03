"""Reader for shapely geometries."""

from __future__ import annotations

from collections.abc import Mapping

from connflex_core.layout import ROLES, Layout, PolygonInput, normalize


def from_shapely(source: Mapping[str, PolygonInput], *, unit: float) -> Layout:
    """Build a layout from a mapping of role name to polygons drawn in ``unit`` meters."""
    unknown = sorted(set(source) - set(ROLES))
    if unknown:
        raise ValueError(f"unknown layout role(s): {', '.join(unknown)}")
    for role in ("device", "handle"):
        if role not in source:
            raise ValueError(f"{role} polygons are required")
    optional = {
        role: normalize(role, source[role], unit) for role in ("release", "vias") if role in source
    }
    return Layout(
        device=normalize("device", source["device"], unit),
        handle=normalize("handle", source["handle"], unit),
        **optional,
    )

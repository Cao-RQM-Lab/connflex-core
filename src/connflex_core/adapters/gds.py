"""Reader for GDSII files through gdstk."""

from __future__ import annotations

import math
from collections.abc import Sequence
from os import PathLike
from typing import TYPE_CHECKING

import shapely
from shapely import GeometryCollection, MultiPolygon, Polygon
from shapely.geometry.base import BaseGeometry

from connflex_core.adapters.shapely_ import from_shapely
from connflex_core.layout import ROLES, LayerMap, Layout, PolygonInput

if TYPE_CHECKING:
    import gdstk


def load_gds(path: str | PathLike[str], *, layers: LayerMap, cell: str | None = None) -> Layout:
    """Read a GDSII file into a layout.

    File units are applied. ``cell`` names the cell to read and is required when a file has
    more than one top-level cell. Polygons with cut lines, which is how GDSII stores holes,
    are repaired into polygons with holes, and each repaired area is verified against the
    file's polygon.
    """
    try:
        import gdstk
    except ImportError as exc:
        raise ImportError(
            "reading GDS files needs gdstk: pip install 'connflex-core[gds]'"
        ) from exc
    library = gdstk.read_gds(str(path))
    top = _top_cell(library, cell, str(path))
    roles: dict[str, PolygonInput] = {}
    for role in ROLES:
        spec: tuple[int, int] | None = getattr(layers, role)
        if spec is None:
            continue
        found = top.get_polygons(
            apply_repetitions=True, include_paths=True, depth=None, layer=spec[0], datatype=spec[1]
        )
        polygons: list[Polygon] = []
        for index, polygon in enumerate(found):
            parts = _polygons(_repair(polygon.points))
            area, expected = sum(part.area for part in parts), polygon.area()
            if not math.isclose(area, expected, rel_tol=1e-6):
                raise ValueError(
                    f"polygon {index} on {role} layer {spec} in cell {top.name!r} of {path} "
                    f"has area {area:g} after repair but {expected:g} in the file"
                )
            polygons.extend(parts)
        if polygons:
            roles[role] = polygons
        elif role in ("device", "handle"):
            raise ValueError(f"cell {top.name!r} of {path} has no polygons on {role} layer {spec}")
    return from_shapely(roles, unit=library.unit)


def _top_cell(library: gdstk.Library, name: str | None, path: str) -> gdstk.Cell:
    import gdstk

    if name is None:
        tops = [top for top in library.top_level() if isinstance(top, gdstk.Cell)]
        if len(tops) != 1:
            names = ", ".join(top.name for top in tops)
            raise ValueError(f"{path} has {len(tops)} top-level cells ({names}); pass cell=")
        return tops[0]
    for candidate in library.cells:
        if isinstance(candidate, gdstk.Cell) and candidate.name == name:
            return candidate
    raise ValueError(f"{path} has no cell named {name!r}")


def _repair(points: Sequence[tuple[float, float]]) -> BaseGeometry:
    return shapely.make_valid(Polygon(points), method="structure", keep_collapsed=False)


def _polygons(geometry: BaseGeometry) -> list[Polygon]:
    if isinstance(geometry, Polygon):
        return [] if geometry.is_empty else [geometry]
    if isinstance(geometry, MultiPolygon | GeometryCollection):
        return [polygon for part in geometry.geoms for polygon in _polygons(part)]
    return []

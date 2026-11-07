"""Reader for GDSII files through gdstk."""

from __future__ import annotations

from os import PathLike
from typing import TYPE_CHECKING

import numpy as np
import shapely
from numpy.typing import NDArray
from shapely import Polygon

from connflex_core.adapters.shapely_ import from_shapely
from connflex_core.geometry.heal import polygon_parts
from connflex_core.layout import ROLES, LayerMap, Layout, PolygonInput

if TYPE_CHECKING:
    import gdstk


def load_gds(path: str | PathLike[str], *, layers: LayerMap, cell: str | None = None) -> Layout:
    """Read a GDSII file into a layout.

    File units are applied. ``cell`` names the cell to read and is required when a file has
    more than one top-level cell. Polygons with cut lines, which is how GDSII stores holes,
    are repaired into polygons with holes, and each repaired area is verified against the
    file's polygon. Polygons are handled as arrays; only invalid ones are rebuilt.
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
        if found:
            shapes = shapely.polygons([shapely.linearrings(polygon.points) for polygon in found])
            repaired = _repair(shapes)
            areas = shapely.area(repaired)
            expected = np.array([polygon.area() for polygon in found])
            for index in np.flatnonzero(~np.isclose(areas, expected, rtol=1e-6, atol=0.0)):
                raise ValueError(
                    f"polygon {index} on {role} layer {spec} in cell {top.name!r} of {path} "
                    f"has area {areas[index]:g} after repair but {expected[index]:g} in the file"
                )
            polygons = [part for shape in repaired for part in polygon_parts(shape)]
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


def _repair(shapes: NDArray[np.object_]) -> NDArray[np.object_]:
    """Rebuild invalid polygons, which is how cut-line holes arrive, and keep valid ones as is."""
    invalid = ~shapely.is_valid(shapes)
    if invalid.any():
        shapes = shapes.copy()
        shapes[invalid] = shapely.make_valid(
            shapes[invalid], method="structure", keep_collapsed=False
        )
    return shapes

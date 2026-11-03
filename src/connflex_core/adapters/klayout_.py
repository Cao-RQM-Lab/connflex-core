"""Reader for KLayout layouts."""

from __future__ import annotations

from typing import TYPE_CHECKING

from shapely import Polygon

from connflex_core.adapters.shapely_ import from_shapely
from connflex_core.layout import ROLES, LayerMap, Layout, PolygonInput

if TYPE_CHECKING:
    import klayout.db as kdb


def from_klayout(
    source: kdb.Layout | kdb.Cell, *, layers: LayerMap, cell: str | None = None
) -> Layout:
    """Read a KLayout layout or cell into a layout.

    Database units are applied and instances are flattened. For a ``Layout``, ``cell`` names
    the cell to read and is required when there is more than one top cell; a ``Cell`` is read
    as given.
    """
    try:
        import klayout.db as kdb
    except ImportError as exc:
        raise ImportError(
            "reading KLayout layouts needs klayout: pip install 'connflex-core[klayout]'"
        ) from exc
    top = source if isinstance(source, kdb.Cell) else _top_cell(source, cell)
    layout = top.layout()
    roles: dict[str, PolygonInput] = {}
    for role in ROLES:
        spec: tuple[int, int] | None = getattr(layers, role)
        if spec is None:
            continue
        index = layout.find_layer(*spec)
        polygons = [] if index is None else _polygons(kdb.Region(top.begin_shapes_rec(index)))
        if polygons:
            roles[role] = polygons
        elif role in ("device", "handle"):
            raise ValueError(f"cell {top.name!r} has no polygons on {role} layer {spec}")
    return from_shapely(roles, unit=layout.dbu * 1e-6)


def _top_cell(layout: kdb.Layout, name: str | None) -> kdb.Cell:
    if name is None:
        tops = layout.top_cells()
        if len(tops) != 1:
            names = ", ".join(top.name for top in tops)
            raise ValueError(f"layout has {len(tops)} top cells ({names}); pass cell=")
        return tops[0]
    found = layout.cell(name)
    if found is None:
        raise ValueError(f"layout has no cell named {name!r}")
    return found


def _polygons(region: kdb.Region) -> list[Polygon]:
    region.merge()
    polygons: list[Polygon] = []
    for polygon in region.each():
        hull = [(point.x, point.y) for point in polygon.each_point_hull()]
        holes = [
            [(point.x, point.y) for point in polygon.each_point_hole(i)]
            for i in range(polygon.holes())
        ]
        polygons.append(Polygon(hull, holes))
    return polygons

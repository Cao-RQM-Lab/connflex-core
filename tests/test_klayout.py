from typing import Any

import pytest

from connflex_core import LayerMap, from_klayout

kdb = pytest.importorskip("klayout.db")

LAYERS = LayerMap(device=(1, 0), handle=(2, 0))


def plate_with_hole() -> Any:
    plate = kdb.Polygon(
        [kdb.Point(0, 0), kdb.Point(200000, 0), kdb.Point(200000, 50000), kdb.Point(0, 50000)]
    )
    plate.insert_hole(
        [
            kdb.Point(90000, 20000),
            kdb.Point(110000, 20000),
            kdb.Point(110000, 40000),
            kdb.Point(90000, 40000),
        ]
    )
    return plate


def test_instances_are_flattened_and_the_database_unit_applied() -> None:
    layout = kdb.Layout()
    layout.dbu = 0.001
    top, child = layout.create_cell("TOP"), layout.create_cell("CHILD")
    child.shapes(layout.layer(1, 0)).insert(plate_with_hole())
    top.insert(kdb.CellInstArray(child.cell_index(), kdb.Trans(kdb.Vector(1_000_000, 0))))
    top.shapes(layout.layer(2, 0)).insert(kdb.Box(-50000, -50000, 1_250_000, 250000))
    result = from_klayout(layout, layers=LAYERS)
    assert result.device.area == pytest.approx((200 * 50 - 20 * 20) * 1e-12)
    assert len(result.device.geoms[0].interiors) == 1
    assert result.device.bounds == pytest.approx((1000e-6, 0, 1200e-6, 50e-6))
    assert result.handle.area == pytest.approx(1300e-6 * 300e-6)
    assert from_klayout(top, layers=LAYERS).device.area == result.device.area


def test_several_top_cells_need_a_name() -> None:
    layout = kdb.Layout()
    for name, size in (("A", 10000), ("B", 20000)):
        cell = layout.create_cell(name)
        cell.shapes(layout.layer(1, 0)).insert(kdb.Box(0, 0, size, size))
        cell.shapes(layout.layer(2, 0)).insert(kdb.Box(0, 0, size, size))
    with pytest.raises(ValueError, match="cell"):
        from_klayout(layout, layers=LAYERS)
    assert from_klayout(layout, layers=LAYERS, cell="B").device.area == pytest.approx(400e-12)

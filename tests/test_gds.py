from pathlib import Path

import pytest

from connflex_core import LayerMap, load_gds

gdstk = pytest.importorskip("gdstk")

LAYERS = LayerMap(device=(1, 0), handle=(2, 0))


@pytest.fixture
def gds_file(tmp_path: Path) -> Path:
    library = gdstk.Library(unit=1e-6, precision=1e-9)
    top = library.new_cell("TOP")
    plate = gdstk.rectangle((0, 0), (200, 50))
    hole = gdstk.rectangle((90, 20), (110, 40))
    top.add(*gdstk.boolean(plate, hole, "not", layer=1, datatype=0))
    top.add(gdstk.rectangle((-50, -50), (250, 250), layer=2, datatype=0))
    path = tmp_path / "layout.gds"
    library.write_gds(path)
    return path


def test_round_trips_areas_and_holes_in_meters(gds_file: Path) -> None:
    layout = load_gds(gds_file, layers=LAYERS)
    assert layout.device.area == pytest.approx((200 * 50 - 20 * 20) * 1e-12)
    assert len(layout.device.geoms) == 1
    assert len(layout.device.geoms[0].interiors) == 1
    assert layout.handle.area == pytest.approx(300 * 300 * 1e-12)
    assert layout.release is None


def test_missing_required_layer_is_an_error(gds_file: Path) -> None:
    with pytest.raises(ValueError, match="handle"):
        load_gds(gds_file, layers=LayerMap(device=(1, 0), handle=(9, 0)))


def test_several_top_cells_need_a_name(tmp_path: Path) -> None:
    library = gdstk.Library(unit=1e-6, precision=1e-9)
    for name, size in (("A", 10), ("B", 20)):
        cell = library.new_cell(name)
        cell.add(gdstk.rectangle((0, 0), (size, size), layer=1, datatype=0))
        cell.add(gdstk.rectangle((0, 0), (size, size), layer=2, datatype=0))
    path = tmp_path / "two.gds"
    library.write_gds(path)
    with pytest.raises(ValueError, match="cell"):
        load_gds(path, layers=LAYERS)
    assert load_gds(path, layers=LAYERS, cell="B").device.area == pytest.approx(400e-12)

from pathlib import Path

import pytest
from shapely import Point, Polygon

from connflex_core import LayerMap, load_gds
from connflex_core.adapters import gds

gdstk = pytest.importorskip("gdstk")

LAYERS = LayerMap(device=(1, 0), handle=(2, 0))
# A truss in micrometers whose open rectangle x -368..-232, y 75..225 is a hole drawn
# with a cut line along y = 75, as GDSII stores holes.
CUT_LINE = [
    (-168, 175), (-150, 175), (-150, 215), (-190, 215), (-190, 175), (-172, 175),
    (-172, 75), (-228, 75), (-228, 225), (-200, 225), (-200, 231), (-400, 231),
    (-400, 225), (-372, 225), (-372, 75), (-428, 75), (-428, 175), (-410, 175),
    (-410, 215), (-450, 215), (-450, 175), (-432, 175), (-432, 75), (-368, 75),
    (-368, 75), (-368, 225), (-232, 225), (-232, 75), (-368, 75), (-368, 75), (-168, 75),
]  # fmt: skip


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


def test_cut_line_along_an_edge_is_repaired_without_filling_the_hole() -> None:
    repaired = gds._repair(CUT_LINE)
    assert repaired.is_valid
    assert repaired.area == pytest.approx(6400, rel=1e-9)
    assert not repaired.covers(Point(-300, 150))


def test_cut_line_hole_round_trips_through_a_file(tmp_path: Path) -> None:
    library = gdstk.Library(unit=1e-6, precision=1e-9)
    top = library.new_cell("TOP")
    top.add(gdstk.Polygon(CUT_LINE, layer=1, datatype=0))
    top.add(gdstk.rectangle((-500, 0), (-100, 300), layer=2, datatype=0))
    path = tmp_path / "cut.gds"
    library.write_gds(path)
    layout = load_gds(path, layers=LAYERS)
    assert layout.device.area == pytest.approx(6400e-12, rel=1e-9)
    assert not layout.device.covers(Point(-300e-6, 150e-6))


def test_repaired_area_is_checked_against_the_file(
    gds_file: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(gds, "_repair", lambda points: Polygon(points).envelope)
    with pytest.raises(ValueError, match=r"polygon 0 on device layer \(1, 0\)"):
        load_gds(gds_file, layers=LAYERS)

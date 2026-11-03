import pytest

from connflex_core import LayerMap, from_gdsfactory

gf = pytest.importorskip("gdsfactory")
kdb = pytest.importorskip("klayout.db")


def test_component_polygons_round_trip_in_meters() -> None:
    gf.gpdk.PDK.activate()
    component = gf.Component()
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
    component.add_polygon(plate, layer=(1, 0))
    component.add_polygon([(-50, -50), (250, -50), (250, 250), (-50, 250)], layer=(2, 0))
    layout = from_gdsfactory(component, layers=LayerMap(device=(1, 0), handle=(2, 0)))
    assert layout.device.area == pytest.approx((200 * 50 - 20 * 20) * 1e-12)
    assert len(layout.device.geoms[0].interiors) == 1
    assert layout.handle.area == pytest.approx(300e-6 * 300e-6)

import pytest
from shapely import Polygon, box

from connflex_core import Layout, adapters


def test_from_shapely_scales_to_meters() -> None:
    layout = Layout.from_shapely(
        device=[box(0, 0, 200, 50)], handle=[box(-50, -50, 250, 250)], unit=1e-6
    )
    assert layout.device.area == pytest.approx(200e-6 * 50e-6)
    assert layout.handle.bounds == pytest.approx((-50e-6, -50e-6, 250e-6, 250e-6))
    assert layout.release is None
    assert layout.vias is None


def test_touching_polygons_are_unioned() -> None:
    layout = Layout.from_shapely(
        device=[box(0, 0, 1, 1), box(1, 0, 2, 1)], handle=box(0, 0, 2, 1), unit=1.0
    )
    assert len(layout.device.geoms) == 1
    assert layout.device.area == pytest.approx(2.0)


def test_invalid_polygon_is_rejected_with_role_and_index() -> None:
    bowtie = Polygon([(0, 0), (1, 1), (1, 0), (0, 1)])
    with pytest.raises(ValueError, match="handle polygon 1"):
        Layout.from_shapely(device=box(0, 0, 1, 1), handle=[box(0, 0, 1, 1), bowtie], unit=1.0)


def test_empty_role_is_rejected() -> None:
    with pytest.raises(ValueError, match="device"):
        Layout.from_shapely(device=[], handle=box(0, 0, 1, 1), unit=1.0)


def test_registry_lists_adapters_and_rejects_unknown_names() -> None:
    assert "shapely" in adapters.available()
    with pytest.raises(KeyError, match="shapely"):
        adapters.get("nope")

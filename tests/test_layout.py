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
    assert {"gds", "gdsfactory", "klayout", "shapely"} <= set(adapters.available())
    with pytest.raises(KeyError, match="shapely"):
        adapters.get("nope")


def test_separate_polygons_keep_their_order_and_only_touching_ones_merge() -> None:
    squares = [box(3 * i, 0, 3 * i + 1, 1) for i in range(50)]
    overlapping = [box(0, 10, 2, 11), box(1, 10, 3, 11)]
    layout = Layout.from_shapely(
        device=[*squares, *overlapping], handle=box(0, 0, 200, 20), unit=1.0
    )
    parts = list(layout.device.geoms)
    assert len(parts) == 51
    assert [part.bounds[0] for part in parts[:50]] == [3 * i for i in range(50)]
    assert parts[50].area == pytest.approx(3.0)


def test_polygons_are_stored_in_canonical_vertex_order() -> None:
    shell = [(0, 0), (10, 0), (10, 10), (0, 10)]
    holes = [[(1, 1), (2, 1), (2, 2), (1, 2)], [(5, 5), (7, 5), (7, 7), (5, 7)]]
    drawn = Polygon(shell, holes)
    redrawn = Polygon(shell[2:] + shell[:2], [holes[1][::-1], holes[0][1:] + holes[0][:1]])
    handle = box(-1, -1, 11, 11)
    first = Layout.from_shapely(device=[drawn], handle=handle, unit=1.0)
    second = Layout.from_shapely(device=[redrawn], handle=handle, unit=1.0)
    assert first.device.equals_exact(second.device, 0)

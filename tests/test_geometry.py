import pytest
from shapely import MultiPolygon, Polygon, box, unary_union

from connflex_core.geometry.bodies import layered_components, select_anchor
from connflex_core.geometry.heal import as_multipolygon, heal, opening
from connflex_core.geometry.release import fill_release_holes, surviving_box

UM = 1e-6


def test_heal_closes_hairline_gaps_and_keeps_real_ones() -> None:
    left = box(0, 0, 10 * UM, 10 * UM)
    near = box(10 * UM + 5e-9, 0, 20 * UM, 10 * UM)
    far = box(20 * UM + 50e-9, 0, 30 * UM, 10 * UM)
    healed = heal(MultiPolygon([left, near, far]), distance=10e-9, tolerance=10e-9)
    assert len(healed.geoms) == 2


def test_opening_removes_narrow_features_and_restores_wide_ones_exactly() -> None:
    pad = box(0, 0, 100 * UM, 100 * UM)
    strip = box(100 * UM, 45 * UM, 300 * UM, 55 * UM)
    opened = opening(as_multipolygon(unary_union([pad, strip])), radius=10 * UM)
    assert opened.area == pytest.approx(pad.area)
    assert opened.bounds == pytest.approx(pad.bounds)


def test_small_holes_are_filled_and_reported() -> None:
    small = box(10 * UM, 10 * UM, 15 * UM, 15 * UM)
    large = box(40 * UM, 40 * UM, 60 * UM, 60 * UM)
    plate = Polygon(box(0, 0, 100 * UM, 100 * UM).exterior, [small.exterior, large.exterior])
    filled, holes = fill_release_holes(MultiPolygon([plate]), max_area=100 * UM * UM)
    assert len(filled.geoms[0].interiors) == 1
    assert [hole.area for hole in holes] == pytest.approx([25 * UM * UM])


def test_surviving_box_excludes_openings_and_undercut() -> None:
    plate = MultiPolygon([box(0, 0, 100 * UM, 100 * UM)])
    opening_ = MultiPolygon([box(40 * UM, 40 * UM, 60 * UM, 60 * UM)])
    survivor = surviving_box(plate, plate, opening_, distance=5 * UM, resolution=8)
    assert survivor.bounds == pytest.approx((5 * UM, 5 * UM, 95 * UM, 95 * UM))
    assert not survivor.intersects(Polygon(box(36 * UM, 50 * UM, 64 * UM, 50.1 * UM)))
    assert survivor.contains(Polygon(box(20 * UM, 48 * UM, 30 * UM, 52 * UM)))


def test_components_connect_only_through_adjacent_layers() -> None:
    top = [box(0, 0, 1, 1), box(5, 5, 6, 6)]
    middle = [box(0.2, 0.2, 0.8, 0.8)]
    bottom = [box(-1, -1, 2, 2), box(4, 4, 7, 7)]
    components = layered_components([top, middle, bottom], seed=(2, 0))
    assert set(components[0]) == {(2, 0), (1, 0), (0, 0)}
    assert sorted(map(sorted, components[1:])) == [[(0, 1)], [(2, 1)]]


def test_anchor_defaults_to_the_largest_handle_polygon_and_accepts_overrides() -> None:
    parts = [box(0, 0, 1, 1), box(10, 10, 20, 20)]
    assert select_anchor(parts, None) == 1
    assert select_anchor(parts, 0) == 0
    assert select_anchor(parts, (0.5, 0.5)) == 0
    with pytest.raises(ValueError, match="anchor"):
        select_anchor(parts, (50.0, 50.0))
    with pytest.raises(ValueError, match="anchor"):
        select_anchor(parts, 2)


def test_opening_treats_each_polygon_on_its_own() -> None:
    def dumbbell(y: float) -> Polygon:
        parts = [box(0, y, 10, y + 10), box(10, y + 4, 20, y + 6), box(20, y, 30, y + 10)]
        merged = unary_union(parts)
        assert isinstance(merged, Polygon)
        return merged

    opened = opening(MultiPolygon([dumbbell(0), dumbbell(100)]), radius=2)
    assert len(opened.geoms) == 4
    assert sorted(part.area for part in opened.geoms) == pytest.approx([100.0] * 4)
    assert opened.bounds == pytest.approx((0, 0, 30, 110))

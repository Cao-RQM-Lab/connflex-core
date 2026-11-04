from collections.abc import Callable

import pytest
from shapely import Point, Polygon, box

from connflex_core import ExtractionParams, Layout, Stack, extract

UM = 1e-6
PAD_A = box(0, 0, 100, 100)
PAD_B = box(300, 0, 400, 100)
BEAM = box(100, 45, 300, 55)
HANDLE_A = box(-50, -50, 150, 150)
Params = Callable[..., ExtractionParams]
LayoutUm = Callable[..., Layout]


def at(x: float, y: float) -> Point:
    return Point(x * UM, y * UM)


def test_clamped_beam_between_an_anchored_and_a_released_pad(
    stack: Stack, params: Params, layout_um: LayoutUm
) -> None:
    structure = extract(layout_um([PAD_A, BEAM, PAD_B], [HANDLE_A]), stack, params())
    anchor, free = structure.bodies
    assert anchor.anchor
    assert anchor.device.contains(at(50, 50))
    assert not anchor.box.is_empty
    assert not free.anchor
    assert free.device.contains(at(350, 50))
    assert free.box.is_empty
    assert free.inertial.mass == pytest.approx(2330 * 25 * UM * (100 * UM) ** 2)
    (beam,) = structure.beams
    assert beam.bodies == (0, 1)
    assert beam.length == pytest.approx(200 * UM, rel=1e-3)
    assert beam.width == pytest.approx(10 * UM, rel=1e-3)
    assert beam.thickness == 25 * UM
    assert beam.start[2] == beam.end[2] == -12.5 * UM


def test_hairline_gap_is_healed(stack: Stack, params: Params, layout_um: LayoutUm) -> None:
    halves = [box(0, 0, 50, 100), box(50.005, 0, 100, 100)]
    structure = extract(layout_um(halves, [HANDLE_A]), stack, params())
    assert len(structure.bodies) == 1
    assert len(structure.bodies[0].device.geoms) == 1


def test_small_holes_become_release_openings_and_large_holes_stay(
    stack: Stack, params: Params, layout_um: LayoutUm
) -> None:
    small = box(10, 10, 15, 15)
    large = box(40, 40, 60, 60)
    plate = Polygon(box(0, 0, 100, 100).exterior, [small.exterior, large.exterior])
    (body,) = extract(layout_um([plate], [HANDLE_A]), stack, params()).bodies
    assert len(body.device.geoms[0].interiors) == 1
    assert body.device.contains(at(12.5, 12.5))
    assert not body.box.intersects(at(12.5, 12.5))
    assert body.box.contains(at(30, 30))


def test_islands_without_a_connection_to_the_anchor_are_dropped(
    stack: Stack, params: Params, layout_um: LayoutUm
) -> None:
    island, island_handle = box(600, 600, 700, 700), box(575, 575, 725, 725)
    structure = extract(layout_um([PAD_A, island], [HANDLE_A, island_handle]), stack, params())
    assert len(structure.bodies) == 1


def test_candidates_touching_one_or_three_bodies_are_rejected(
    stack: Stack, params: Params, layout_um: LayoutUm
) -> None:
    stub = box(100, 45, 150, 55)
    structure = extract(layout_um([PAD_A, stub], [HANDLE_A]), stack, params())
    assert structure.beams == ()
    (candidate,) = structure.rejected
    assert "1 bod" in candidate.reason
    assert candidate.bodies == (0,)

    pad_c = box(150, 300, 250, 400)
    arms = [box(100, 45, 300, 55), box(195, 55, 205, 300)]
    structure = extract(layout_um([PAD_A, PAD_B, pad_c, *arms], [HANDLE_A]), stack, params())
    assert len(structure.bodies) == 3
    assert structure.beams == ()
    junctions = [candidate for candidate in structure.rejected if "3 bod" in candidate.reason]
    assert [candidate.bodies for candidate in junctions] == [(0, 1, 2)]
    assert [0, 1, 2] in [entry["bodies"] for entry in structure.to_dict()["rejected"]]


def test_beam_width_bounds_reject_slivers(
    stack: Stack, params: Params, layout_um: LayoutUm
) -> None:
    thin = box(100, 49, 300, 51)
    assert len(extract(layout_um([PAD_A, thin, PAD_B], [HANDLE_A]), stack, params()).beams) == 1
    structure = extract(
        layout_um([PAD_A, thin, PAD_B], [HANDLE_A]), stack, params(beam_width_min=5 * UM)
    )
    assert structure.beams == ()
    assert any("beam_width_min" in candidate.reason for candidate in structure.rejected)


def test_anchor_can_be_chosen_by_point(stack: Stack, params: Params, layout_um: LayoutUm) -> None:
    handle_b = box(250, -50, 450, 150)
    structure = extract(
        layout_um([PAD_A, BEAM, PAD_B], [HANDLE_A, handle_b]),
        stack,
        params(anchor=(350 * UM, 50 * UM)),
    )
    assert structure.anchor.device.contains(at(350, 50))
    assert len(structure.bodies) == 2
    assert len(structure.beams) == 1


def test_release_openings_are_removed_from_the_device_mass(
    stack: Stack, params: Params, layout_um: LayoutUm
) -> None:
    hole = box(340, 40, 345, 45)
    solid_mass = 2330 * 25 * UM * (100 * UM) ** 2
    hole_mass = 2330 * 25 * UM * (5 * UM) ** 2

    perforated = Polygon(PAD_B.exterior, [hole.exterior])
    structure = extract(layout_um([PAD_A, BEAM, perforated], [HANDLE_A]), stack, params())
    assert structure.bodies[1].inertial.mass == pytest.approx(solid_mass - hole_mass)

    structure = extract(
        layout_um([PAD_A, BEAM, PAD_B], [HANDLE_A], release=[hole]), stack, params()
    )
    assert structure.bodies[1].inertial.mass == pytest.approx(solid_mass - hole_mass)

    straddling = box(390, 40, 410, 45)
    structure = extract(
        layout_um([PAD_A, BEAM, PAD_B], [HANDLE_A], release=[straddling]), stack, params()
    )
    inside_mass = 2330 * 25 * UM * (10 * UM) * (5 * UM)
    assert structure.bodies[1].inertial.mass == pytest.approx(solid_mass - inside_mass)

    corner = box(380, 80, 388, 88)
    perforated = Polygon(PAD_B.exterior, [corner.exterior])
    structure = extract(layout_um([PAD_A, BEAM, perforated], [HANDLE_A]), stack, params())
    pad_area, hole_area = 100 * 100, 8 * 8
    expected_x = (350 * pad_area - 384 * hole_area) / (pad_area - hole_area) * UM
    expected_y = (50 * pad_area - 84 * hole_area) / (pad_area - hole_area) * UM
    assert structure.bodies[1].inertial.center[:2] == pytest.approx([expected_x, expected_y])


def test_candidates_joined_by_a_sliver_narrower_than_heal_distance_are_separate_beams(
    stack: Stack, params: Params, layout_um: LayoutUm
) -> None:
    legs = [box(100, 15, 300, 25), box(100, 75, 300, 85)]
    bridge = box(200, 25, 200.002, 75)
    structure = extract(layout_um([PAD_A, *legs, bridge, PAD_B], [HANDLE_A]), stack, params())
    assert len(structure.beams) == 2
    assert all(beam.width == pytest.approx(10 * UM, rel=1e-3) for beam in structure.beams)
    assert all(beam.length == pytest.approx(200 * UM, rel=1e-3) for beam in structure.beams)


def test_cavity_that_healing_closes_is_a_release_hole(
    stack: Stack, params: Params, layout_um: LayoutUm
) -> None:
    cavity, slit = box(40, 40, 45, 45), box(45, 42, 100, 42.015)
    pad = PAD_A.difference(cavity.union(slit))
    structure = extract(layout_um([pad, BEAM, PAD_B], [HANDLE_A]), stack, params())
    anchor = structure.bodies[0]
    assert len(anchor.device.geoms[0].interiors) == 0
    assert not anchor.box.intersects(at(42.5, 42.5))
    assert anchor.box.contains(at(30, 30))
    device_mass = anchor.inertial.mass - 2330 * 400 * UM * anchor.handle.area
    device_mass -= 2200 * 2 * UM * anchor.box.area
    assert device_mass == pytest.approx(2330 * 25 * UM * (100 * 100 - 5 * 5) * UM * UM, rel=1e-4)

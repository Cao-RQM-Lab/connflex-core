from dataclasses import replace

import numpy as np
import pytest
from shapely import box

from connflex_core import (
    ElectrostaticParams,
    ExtractionParams,
    Layer,
    Layout,
    Material,
    Stack,
    Structure,
    extract,
)
from connflex_core.electrostatics import fringe
from connflex_core.electrostatics.section import Cut, solve_section

UM = 1e-6
BOX = 2 * UM
PITCH, WIDTH, GAP, FINGERS = 14.0, 4.0, 3.0, 4


def comb_over_handle(thickness: float = 25.0) -> Structure:
    """Bonded pad with five fingers and a released shuttle with four between them, on
    two beams, over the handle, gaps of 3 um."""
    stator = [box(300, 100 + PITCH * k, 360, 100 + PITCH * k + WIDTH) for k in range(FINGERS + 1)]
    rotor = [
        box(320, 100 + PITCH * k + WIDTH + GAP, 380, 100 + PITCH * k + 2 * WIDTH + GAP)
        for k in range(FINGERS)
    ]
    spine, beams = box(380, 60, 420, 200), [box(420, 70, 620, 80), box(420, 180, 620, 190)]
    device = [box(0, 0, 300, 300), *stator, *rotor, spine, *beams, box(620, 0, 920, 300)]
    layout = Layout.from_shapely(device=device, handle=[box(-50, -50, 970, 350)], unit=UM)
    stack = Stack(
        device=Layer(thickness=thickness * UM, material=Material.silicon()),
        box=Layer(thickness=BOX, material=Material.oxide()),
        handle=Layer(thickness=400 * UM, material=Material.silicon()),
    )
    params = ExtractionParams(
        release_hole_max_area=100 * UM * UM, release_distance=60 * UM, body_min_width=20 * UM
    )
    return extract(layout, stack, params)


FIELD = ElectrostaticParams(
    element_size=4 * UM,
    interaction_distance=30 * UM,
    dense_cutoff=float("inf"),
    plates_cutoff=float("inf"),
)


@pytest.fixture(scope="module")
def comb() -> Structure:
    return comb_over_handle()


def electrodes(structure: Structure) -> dict[str, int]:
    return {
        "stator": structure.electrode_at((150 * UM, 150 * UM), "device").index,
        "shuttle": structure.electrode_at((400 * UM, 130 * UM), "device").index,
        "handle": structure.electrode_at((400 * UM, 20 * UM), "handle").index,
    }


def walls(handle: bool, offset: float = 0.0) -> Cut:
    """Two wide blocks across a gap of 4 um, the second one ``offset`` farther."""
    start = np.array([-30.0, 4.0 + offset]) * UM
    end = np.array([0.0, 26.0 + offset]) * UM
    body = np.array([0, 1, 2] if handle else [0, 1])
    return Cut(start, end, np.array([0, 1]), body, handle, True)


def test_cross_section_excess_stays_bounded_as_walls_grow_tall() -> None:
    """The fringe at the top and bottom of a wall does not grow with its height, so it
    matters less and less against the field across the gap."""
    excess = []
    for thickness in (10.0, 20.0, 40.0, 80.0):
        section = solve_section(walls(False), thickness * UM, BOX, 30 * UM, 16.0)
        excess.append(-section.value[1])
        assert section.value.sum() == pytest.approx(0.0, abs=1e-9)
    assert all(e > 0 for e in excess)
    assert excess[3] == pytest.approx(excess[2], rel=0.05)
    across = np.array([10.0, 20.0, 40.0, 80.0]) / 4.0
    share = np.array(excess) / across
    assert share[3] < 0.55 * share[2]


def test_cross_section_derivatives_match_finite_differences() -> None:
    """Travel of the far wall along the normal and rise of the near one, with the
    handle beneath; without it the field is even about mid-height and pulls no way."""
    level = solve_section(walls(False), 25 * UM, BOX, 30 * UM, 16.0)
    assert level.gradient[:, 0, 1] == pytest.approx(0.0, abs=1e-9 * np.max(abs(level.gradient)))
    assert level.hessian[1, 0, 1, 0, 1] > 0
    base = solve_section(walls(True), 25 * UM, BOX, 30 * UM, 16.0)
    step = 1e-3
    ahead = solve_section(walls(True, step), 25 * UM, BOX, 30 * UM, 16.0)
    behind = solve_section(walls(True, -step), 25 * UM, BOX, 30 * UM, 16.0)
    step *= UM
    slope = (ahead.value - behind.value) / (2 * step)
    assert base.gradient[:, 1, 0] == pytest.approx(slope, rel=1e-6, abs=1e-9 * np.max(abs(slope)))
    bend = (ahead.gradient - behind.gradient) / (2 * step)
    assert base.hessian[:, :, :, 1, 0] == pytest.approx(
        bend, rel=1e-5, abs=1e-8 * np.max(abs(bend))
    )


def test_fringe_is_a_symmetric_part_of_the_matrix(comb: Structure) -> None:
    model = comb.capacitance(FIELD)
    at = electrodes(comb)
    part = model.fringe
    assert part == pytest.approx(part.T)
    assert np.max(abs(part.sum(axis=1))) < 1e-12 * np.max(abs(part))
    assert model.matrix == pytest.approx(model.oxide + model.released + model.sidewall + part)
    assert -part[at["stator"], at["shuttle"]] > 0
    assert -part[at["shuttle"], at["handle"]] > 0
    assert model.clusters[0].sections > 0


def test_fringe_derivatives_are_slopes_of_the_model(comb: Structure) -> None:
    """At rest and with the shuttle moved, tilted, and turned. Walls that a cut meets
    behind its own, of a third body, move as if square to it, which leaves second
    derivatives in turn a little off."""
    at = electrodes(comb)
    rest = comb.capacitance(FIELD)
    potentials = np.zeros(len(comb.electrodes))
    potentials[at["stator"]], potentials[at["handle"]] = 10.0, 5.0
    for place in (np.zeros(6), np.array([2e-6, 0.4e-6, -0.3e-6, 2e-4, -1e-4, 3e-4])):
        here = rest.at(place)
        stiffness = here.stiffness(potentials)
        for dof in range(6):
            step = np.zeros(6)
            step[dof] = 1e-10 if dof < 3 else 1e-8
            ahead, behind = rest.at(place + step), rest.at(place - step)
            slope = (ahead.matrix - behind.matrix) / (2 * step[dof])
            expected = here.gradient[:, :, dof]
            floor = 1e-6 * np.max(abs(slope)) + 1e-12 * np.max(abs(here.gradient))
            assert expected == pytest.approx(slope, rel=1e-5, abs=floor)
            softening = -(ahead.force(potentials) - behind.force(potentials)) / (2 * step[dof])
            scale = np.sqrt(abs(np.diag(stiffness) * stiffness[dof, dof]))
            assert np.all(abs(stiffness[:, dof] - softening) < 2e-3 * abs(softening) + 1e-5 * scale)


def test_faced_body_is_carried_by_the_gap_ahead_of_the_wall() -> None:
    """A cut across a wall facing a finger whose own silicon lies nearer behind the wall:
    the faced body is expanded from the end of its silicon ahead, whose rate the
    expansion applies, not from its nearest end, equation (160)."""
    plan = fringe.FringePlan((), np.zeros(0, dtype=np.intp), 25 * UM, BOX, 30 * UM, 16.0)
    spans = [[-5.4 * UM, -3.4 * UM, 1], [-2 * UM, 0.0, 0], [4.4 * UM, 6.4 * UM, 1]]
    cutting = fringe._keyed(spans, 0, 1, np.array([0, 1]), -1, plan, 2)
    assert cutting.faced == 1
    assert cutting.position[1] == pytest.approx(4.4 * UM)
    assert cutting.offset[1] == pytest.approx(-0.1 * UM)
    assert cutting.position[0] == 0.0  # the wall itself


def test_handle_lifts_fingers_and_alignment_holds_without_it(comb: Structure) -> None:
    """A grounded shuttle between fingers at 10 V over a grounded handle is pushed up:
    the handle draws the field under the stator fingers down. With no handle beneath,
    the field is even about mid-height, and leaving the stator's height costs field."""
    at = electrodes(comb)
    model = comb.capacitance(FIELD)
    up = model.dofs.index((1, "uz"))
    potentials = np.zeros(len(comb.electrodes))
    potentials[at["stator"]] = 10.0
    assert model.force(potentials)[up] > 0
    alone = comb.capacitance(replace(FIELD, include_plates=False))
    assert alone.fringe[:, at["handle"]] == pytest.approx(0.0)
    assert abs(alone.force(potentials)[up]) < 1e-9 * abs(model.force(potentials)[up])
    assert alone.stiffness(potentials)[up, up] > 0


def test_fringe_share_halves_as_the_device_layer_doubles() -> None:
    """Against the field across the gaps, which grows with thickness, the fringe at the
    top and bottom of the walls stays nearly the same."""
    shares = []
    for thickness in (25.0, 50.0):
        structure = comb_over_handle(thickness)
        at = electrodes(structure)
        model = structure.capacitance(replace(FIELD, include_plates=False))
        pair = (at["stator"], at["shuttle"])
        shares.append(model.fringe[pair] / model.sidewall[pair])
    assert 0 < shares[1] < 0.6 * shares[0]
    assert shares[1] > 0.4 * shares[0]


def test_pull_stays_continuous_as_one_tooth_passes_another() -> None:
    """A tooth on a bar slides under a tooth of an anchored block, 4 um across. Ends of
    facing walls that pass each other round the overlap over the box thickness, so the
    pull from the fringe changes without jumps."""
    stack = Stack(
        device=Layer(thickness=25 * UM, material=Material.silicon()),
        box=Layer(thickness=BOX, material=Material.oxide()),
        handle=Layer(thickness=400 * UM, material=Material.silicon()),
    )
    device = [box(-200, 0, 200, 40), box(-5, 40, 5, 50), box(-150, 64, 150, 264)]
    device += [box(-5, 54, 5, 64), box(-199, -800, -197, 0), box(197, -800, 199, 0)]
    device += [box(-350, -1100, -150, -800), box(150, -1100, 350, -800)]
    layout = Layout.from_shapely(device=device, handle=[box(-400, -1200, 400, 400)], unit=UM)
    extraction = ExtractionParams(
        release_hole_max_area=100 * UM * UM, release_distance=60 * UM, body_min_width=20 * UM
    )
    structure = extract(layout, stack, extraction)
    model = structure.capacitance(
        replace(FIELD, include_plates=False, interaction_distance=10 * UM)
    )
    along = model.dofs.index((model.dofs[0][0], "ux"))
    stator = structure.electrode_at((0.0, 150 * UM), "device").index
    potentials = np.zeros(len(structure.electrodes))
    potentials[stator] = 50.0

    def pull(travel: float) -> float:
        place = np.zeros(len(model.dofs))
        place[along] = travel * UM
        return float(model.at(place, exact=True).force(potentials)[along])

    largest = max(abs(pull(travel)) for travel in (-6.0, -3.0, -1.0, 1.0, 3.0, 6.0))
    for passing in (0.0, 10.0):
        step = 1e-3
        assert abs(pull(passing + step) - pull(passing - step)) < 1e-3 * largest
    unsmoothed = 2 * 0.5 * 50.0**2 * 8.8541878188e-12 * 0.3
    assert unsmoothed > 1e-3 * largest

from dataclasses import replace

import numpy as np
import pytest
from shapely import box

from connflex_core import (
    CapacitanceModel,
    ElectrostaticParams,
    Layout,
    PiecewiseLinear,
    Source,
    Structure,
    TransientOptions,
    extract,
)
from test_equilibrium import EXTRACTION, STACK, UM, net_at, shuttle_layout

FIELD = ElectrostaticParams(
    element_size=8 * UM,
    corner_grading=1.0,
    interaction_distance=10 * UM,
    include_plates=False,
    dense_cutoff=float("inf"),
    plates_cutoff=float("inf"),
)


@pytest.fixture(scope="module")
def shuttle() -> Structure:
    return extract(shuttle_layout(), STACK, EXTRACTION)


@pytest.fixture(scope="module")
def exact(shuttle: Structure) -> CapacitanceModel:
    return shuttle.capacitance(FIELD)


@pytest.fixture(scope="module")
def soft() -> Structure:
    """The comb of the shuttle on two long beams across its fingers, soft along them."""
    device = [box(-20, -110, 20, 110)]
    for side in (-1, 1):
        low, high = sorted((side * 100, side * 400))
        device.append(box(low, -70, high, 70))
        for k in range(5):
            low, high = sorted((side * 40, side * 100))
            device.append(box(low, -30 + 14 * k, high, -26 + 14 * k))
        for k in range(4):
            low, high = sorted((side * 20, side * 80))
            device.append(box(low, -23 + 14 * k, high, -19 + 14 * k))
    for sign in (-1, 1):
        low, high = sorted((sign * 110, sign * 510))
        device.append(box(-2, low, 2, high))
        low, high = sorted((sign * 510, sign * 820))
        device.append(box(-150, low, 150, high))
    layout = Layout.from_shapely(device=device, handle=[box(-600, -900, 600, 900)], unit=UM)
    return extract(layout, STACK, EXTRACTION)


def path() -> list[np.ndarray]:
    """Travel along the fingers and toward one side, where capacitance curves."""
    places = []
    for step in np.linspace(0.0, 1.0, 31):
        place = np.zeros(6)
        place[0], place[1] = 2e-6 * step, 0.6e-6 * step**2
        places.append(place)
    return places


def relative(found: np.ndarray, truth: np.ndarray) -> float:
    return float(np.linalg.norm(found - truth) / np.linalg.norm(truth))


def test_surrogate_follows_field_solution_within_its_tolerance(
    shuttle: Structure, exact: CapacitanceModel
) -> None:
    """Along a path that closes side gaps by a fifth, the expansion of the nearest field
    solution keeps gradient and capacitance within the tolerance, which is a target the
    radius rule aims at, from fewer solutions the looser the tolerance."""
    truths = [exact.at(place) for place in path()]
    anchors = {}
    for tolerance in (1e-2, 3e-4):
        model = shuttle.capacitance(replace(FIELD, surrogate_tolerance=tolerance))
        worst = 0.0
        for place, truth in zip(path(), truths, strict=True):
            found = model.at(place)
            worst = max(
                worst,
                relative(found.gradient, truth.gradient),
                relative(found.matrix, truth.matrix),
            )
            assert found.clearance == pytest.approx(truth.clearance)
        anchors[tolerance] = found.clusters[0].anchors
        assert worst < tolerance
    assert 1 <= anchors[1e-2] < anchors[3e-4] < len(path())
    assert exact.clusters[0].anchors == 0


def test_surrogate_holds_its_tolerance_after_rise_and_tilt(shuttle: Structure) -> None:
    """Rise and tilt change the fringe but not the sidewalls: steps out of the plane must
    not grow a radius that later admits in-plane steps beyond the tolerance."""
    tolerance = 1e-4
    field = replace(FIELD, include_plates=True, interaction_distance=20 * UM)
    exact = shuttle.capacitance(field)
    model = shuttle.capacitance(replace(field, surrogate_tolerance=tolerance))
    places = []
    for rise in (-0.1, -0.2, -0.3, -0.4):  # each within four radii of the last
        place = np.zeros(len(exact.dofs))
        place[2], place[3] = rise * UM, 1e-4
        places.append(place)
    place = places[-1].copy()
    place[0] = 0.2 * UM  # within the radius that the rise steps would have grown
    places.append(place)
    for place in places:
        found, truth = model.at(place), exact.at(place)
        assert relative(found.matrix, truth.matrix) < 3 * tolerance
        assert relative(found.gradient, truth.gradient) < 3 * tolerance
        assert relative(found.fringe, truth.fringe) < 3 * tolerance


def test_solution_at_a_kept_motion_replaces_the_kept_one(shuttle: Structure) -> None:
    model = shuttle.capacitance(replace(FIELD, surrogate_tolerance=1e-2))
    assert model.clusters[0].anchors == 1
    rest = np.zeros(len(model.dofs))
    assert model.at(rest, exact=True).clusters[0].anchors == 1
    place = path()[-1]
    kept = model.at(place).clusters[0].anchors
    assert model.at(place, exact=True).clusters[0].anchors == kept
    assert model.at(place, exact=True).clusters[0].anchors == kept


def test_expanded_solution_at_a_still_motion_is_not_kept_as_the_rest_field(
    shuttle: Structure, exact: CapacitanceModel
) -> None:
    """A surrogate hit at a motion without in-plane travel is an expansion of a kept
    solution; an exact evaluation at rest afterwards must solve the field or take the
    solution found at rest, never that expansion."""
    model = shuttle.capacitance(replace(FIELD, surrogate_tolerance=1e-2))
    first = np.zeros(len(model.dofs))
    first[0], first[2] = 0.2 * UM, 1.0 * UM
    second = np.zeros(len(model.dofs))
    second[2] = 1.0 * UM
    moved = model.at(first).at(second)
    assert moved.clusters[0].anchors == 2  # a hit on the solution at ``first``
    rest = moved.at(np.zeros(len(model.dofs)), exact=True)
    assert rest.sidewall.tolist() == exact.sidewall.tolist()


def test_exact_evaluation_bypasses_surrogate(shuttle: Structure, exact: CapacitanceModel) -> None:
    model = shuttle.capacitance(replace(FIELD, surrogate_tolerance=1e-2))
    place = path()[-1]
    model.at(path()[-2])
    truth, found = exact.at(place), model.at(place, exact=True)
    assert found.gradient == pytest.approx(truth.gradient, rel=1e-12, abs=1e-30)
    assert found.matrix == pytest.approx(truth.matrix, rel=1e-12, abs=1e-30)


def test_equilibrium_with_surrogate_matches_exact_and_reports_exact_stiffness(
    soft: Structure,
) -> None:
    """120 V pulls the comb in by a sixth of its gap; the equilibrium found through the
    surrogate is that of the field solution, and its model is the exact one there."""
    exact = soft.capacitance(FIELD)
    model = soft.capacitance(replace(FIELD, surrogate_tolerance=1e-4))
    potentials = {net_at(soft, 250, 0): 120.0}
    truth = soft.equilibrium(exact, potentials)
    found = soft.equilibrium(model, potentials)
    assert truth.displacement[0] > 0.1 * 3 * UM
    assert found.displacement == pytest.approx(
        truth.displacement, rel=1e-6, abs=1e-6 * truth.displacement[0]
    )
    there = exact.at(found.displacement)
    assert found.model.gradient == pytest.approx(there.gradient, rel=1e-12, abs=1e-30)
    assert found.stiffness == pytest.approx(
        truth.stiffness, rel=1e-5, abs=1e-5 * float(np.max(abs(truth.stiffness)))
    )


def test_transient_with_surrogate_settles_with_few_field_solutions(soft: Structure) -> None:
    model = soft.capacitance(replace(FIELD, surrogate_tolerance=1e-4))
    right = soft.electrode_at((250 * UM, 0.0), "device")
    slowest = 2 * np.pi * soft.modal().frequencies_hz[0]
    ramp = PiecewiseLinear([0.0, 2 / slowest], [0.0, 120.0])
    times = np.linspace(0.0, 60 / slowest, 61)
    run = soft.transient(model, [Source(right, ramp)], times, TransientOptions(quality_factor=2.0))
    rest = soft.equilibrium(soft.capacitance(FIELD), {right.net: 120.0})
    assert run.displacement[-1] == pytest.approx(
        rest.displacement, rel=1e-4, abs=1e-4 * rest.displacement[0]
    )
    assert run.model.clusters[0].anchors < run.evaluations / 20


def test_surrogate_tolerance_is_checked() -> None:
    with pytest.raises(ValueError, match="surrogate_tolerance"):
        ElectrostaticParams(surrogate_tolerance=0.0)
    with pytest.raises(ValueError, match="surrogate_tolerance"):
        ElectrostaticParams(surrogate_tolerance=1.0)


def test_surrogate_hit_reuses_the_mapped_hessian_blocks(shuttle: Structure) -> None:
    """Expansion leaves the second derivative of a kept solution as it is, so its map to
    the degrees of freedom is built once per kept solution and shared by every hit."""
    model = shuttle.capacitance(replace(FIELD, surrogate_tolerance=1e-2))
    place = path()[-2]
    first = model.at(place)
    second = model.at(place * 1.001)
    assert second.clusters[0].anchors == first.clusters[0].anchors
    assert len(second.blocks) == len(first.blocks) > 0
    assert all(
        mine.values is theirs.values
        for mine, theirs in zip(second.blocks, first.blocks, strict=True)
    )


def test_expand_follows_at_to_second_order_with_the_clearance_of_at(
    exact: CapacitanceModel,
) -> None:
    """A model expanded at a nearby displacement, equation (149) over the whole model,
    matches the model solved there to second order in the matrix and first in the
    gradient, so halving the step cuts the errors by about eight and four, and its
    clearance is that of the solved model."""
    matrix_errors, gradient_errors = [], []
    for scale in (1.0, 0.5):
        place = np.zeros(len(exact.dofs))
        place[0], place[1] = 0.3e-6 * scale, 0.15e-6 * scale
        truth, near = exact.at(place), exact.expand(place)
        matrix_errors.append(relative(near.matrix, truth.matrix))
        gradient_errors.append(relative(near.gradient, truth.gradient))
        assert near.displacement.tolist() == place.tolist()
        assert near.clearance == truth.clearance
        assert near.closest == truth.closest
        assert near.blocks is exact.blocks
    assert 6 < matrix_errors[0] / matrix_errors[1] < 10
    assert 3 < gradient_errors[0] / gradient_errors[1] < 5
    same = exact.expand(np.zeros(len(exact.dofs)))
    assert same.matrix.tolist() == exact.matrix.tolist()
    assert same.gradient.tolist() == exact.gradient.tolist()

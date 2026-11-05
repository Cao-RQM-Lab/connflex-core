from dataclasses import replace

import numpy as np
import pytest
import scipy.optimize
from shapely import Polygon, box

from connflex_core import (
    ElectrostaticParams,
    Equilibrium,
    ExtractionParams,
    Layer,
    Layout,
    Material,
    Stack,
    StaticOptions,
    Structure,
    extract,
)

UM = 1e-6
EPS0 = 8.8541878188e-12
THICKNESS, BOX = 25 * UM, 2 * UM
SIDE = 100.0
GAP, WIDTH, OVERLAP, GAPS = 3.0, 4.0, 40.0, 8
STACK = Stack(
    device=Layer(thickness=THICKNESS, material=Material.silicon()),
    box=Layer(thickness=BOX, material=Material.oxide()),
    handle=Layer(thickness=400 * UM, material=Material.silicon()),
)
EXTRACTION = ExtractionParams(
    release_hole_max_area=100 * UM * UM, release_distance=60 * UM, body_min_width=20 * UM
)
HANDLE_ALONE = ElectrostaticParams(include_sidewalls=False)
WALLS_ALONE = ElectrostaticParams(
    sidewall_method="plates", interaction_distance=10 * UM, include_plates=False
)


def plate_on_beams(*more: Polygon, sliver: float = 0.0) -> Structure:
    """Released square plate over the handle on four beams, two to a pad on either side."""
    beams = [box(x, y, x + 200, y + 10) for x in (-200, SIDE) for y in (10, SIDE - 20)]
    pads = [box(-500, -100, -200, 200), box(SIDE + 200, -100, SIDE + 500, 200)]
    layout = Layout.from_shapely(
        device=[box(0, 0, SIDE, SIDE), *beams, *pads, *more],
        handle=[box(-600, -200, 700, 300)],
        unit=UM,
    )
    return extract(layout, STACK, replace(EXTRACTION, beam_width_min=sliver * UM))


@pytest.fixture(scope="module")
def plate() -> Structure:
    return plate_on_beams()


@pytest.fixture(scope="module")
def shuttle() -> Structure:
    """Spine on four beams with a comb of eight gaps toward a stator on either side."""
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
        for y in (-105, 95):
            low, high = sorted((side * 20, side * 220))
            device.append(box(low, y, high, y + 10))
            low, high = sorted((side * 220, side * 520))
            device.append(
                box(low, y - 290 if y < 0 else y - 10, high, y + 20 if y < 0 else y + 300)
            )
    layout = Layout.from_shapely(device=device, handle=[box(-600, -450, 600, 450)], unit=UM)
    return extract(layout, STACK, EXTRACTION)


def net_at(structure: Structure, x: float, y: float, layer: str = "device") -> int:
    return structure.electrode_at((x * UM, y * UM), layer).net


def test_plate_settles_where_its_spring_balances_the_pull(plate: Structure) -> None:
    stiffness = plate.modal().stiffness
    spring, area, volts = stiffness[2, 2], (SIDE * UM) ** 2, 300.0
    assert np.max(abs(np.delete(stiffness[2], 2))) < 1e-9 * spring

    def unbalance(rise: float) -> float:
        return float(spring * rise + 0.5 * EPS0 * area * volts**2 / (BOX + rise) ** 2)

    rise = scipy.optimize.brentq(unbalance, -BOX / 3, 0.0, xtol=1e-18)
    model = plate.capacitance(HANDLE_ALONE)
    handle = net_at(plate, 0, 0, "handle")
    result = plate.equilibrium(model, {handle: volts})
    assert isinstance(result, Equilibrium)
    assert result.dofs == model.dofs
    assert result.displacement[2] == pytest.approx(rise, rel=1e-8)
    assert np.max(abs(np.delete(result.displacement, 2))) < 1e-9 * abs(rise)
    assert result.force[2] == pytest.approx(spring * rise, rel=1e-8)
    softened = spring - EPS0 * area * volts**2 / (BOX + rise) ** 3
    assert result.stiffness[2, 2] == pytest.approx(softened, rel=1e-8)
    assert result.stable
    assert result.softest > 0
    assert result.model.displacement == pytest.approx(result.displacement)
    under = plate.electrode_at((0.0, 0.0), "handle").index
    assert result.potentials[under] == volts
    assert result.charge == pytest.approx(result.model.matrix @ result.potentials)
    assert result.charge[under] > 0


def test_comb_travel_grows_with_voltage_squared(shuttle: Structure) -> None:
    spring = shuttle.modal().stiffness[0, 0]
    model = shuttle.capacitance(WALLS_ALONE)
    right = net_at(shuttle, 250, 0)
    travel = [shuttle.equilibrium(model, {right: volts}).displacement for volts in (50.0, 100.0)]
    pull = GAPS * EPS0 * THICKNESS * 50.0**2 / (2 * GAP * UM)
    assert travel[0][0] == pytest.approx(pull / spring, rel=1e-8)
    assert travel[1][0] == pytest.approx(4 * travel[0][0], rel=1e-8)
    assert np.max(abs(travel[1][1:])) < 1e-9 * travel[1][0]


def test_plate_is_followed_up_to_pull_in_at_the_closed_form_voltage(plate: Structure) -> None:
    spring, area = plate.modal().stiffness[2, 2], (SIDE * UM) ** 2
    pull_in = np.sqrt(8 * spring * BOX**3 / (27 * EPS0 * area))
    model = plate.capacitance(HANDLE_ALONE)
    handle = net_at(plate, 0, 0, "handle")

    below = plate.equilibrium(model, {handle: 0.9 * pull_in}, options=StaticOptions(steps=10))
    assert (below.scale, below.pull_in, below.limit) == (1.0, None, None)
    assert below.ramp_scales.tolist() == [k / 10 for k in range(11)]
    for scale, place in zip(below.ramp_scales, below.ramp_displacements, strict=True):
        volts = scale * 0.9 * pull_in
        balance = spring * place[2] + 0.5 * EPS0 * area * volts**2 / (BOX + place[2]) ** 2
        assert abs(balance) < 1e-8 * spring * BOX
    assert below.ramp_displacements[-1] == pytest.approx(below.displacement)
    assert np.all(np.diff(below.ramp_displacements[:, 2]) < 0)
    assert below.ramp_softest[-1] == below.softest

    beyond = plate.equilibrium(
        model, {handle: 1.2 * pull_in}, options=StaticOptions(pull_in_resolution=1e-6)
    )
    assert beyond.limit in ("unstable", "no equilibrium")
    assert beyond.pull_in is not None
    assert beyond.pull_in * 1.2 == pytest.approx(1.0, rel=1e-5)
    assert 0 < beyond.pull_in - beyond.scale <= 1e-6
    assert beyond.stable
    assert np.max(beyond.potentials) == pytest.approx(beyond.scale * 1.2 * pull_in)
    assert beyond.displacement[2] == pytest.approx(-BOX / 3, rel=1e-2)


def test_comb_loses_sideways_stability_at_the_closed_form_voltage(shuttle: Structure) -> None:
    """Both stators at one potential pull the shuttle equally: it stays centered until
    the spring across the fingers is used up."""
    spring = shuttle.modal().stiffness[1, 1]
    limit = np.sqrt(spring * (GAP * UM) ** 3 / (2 * GAPS * EPS0 * THICKNESS * OVERLAP * UM))
    model = shuttle.capacitance(WALLS_ALONE)
    stators = {net_at(shuttle, x, 0): 1.3 * limit for x in (-250, 250)}
    result = shuttle.equilibrium(model, stators, options=StaticOptions(pull_in_resolution=1e-5))
    assert result.limit == "unstable"
    assert result.pull_in is not None
    assert result.pull_in * 1.3 == pytest.approx(1.0, rel=1e-4)
    assert result.stable
    assert np.max(abs(result.displacement[:3])) < 1e-9 * GAP * UM


def test_mechanical_load_adds_to_electrostatic_force(plate: Structure) -> None:
    model = plate.capacitance(HANDLE_ALONE)
    load = {1: [0.0, 0.0, 1e-6, 0.0, 0.0, 0.0]}
    result = plate.equilibrium(model, {}, loads=load)
    assert result.displacement == pytest.approx(plate.static(load), rel=1e-9, abs=1e-20)
    assert not result.force.any()


def test_body_that_no_beam_holds_stays_at_rest(plate: Structure) -> None:
    """A released island hangs from a pad by a neck too thin to count as a beam."""
    islanded = plate_on_beams(box(-380, 240, -340, 280), box(-361, 200, -359, 240), sliver=5.0)
    assert len(islanded.bodies) == 3
    assert len(islanded.beams) == 4
    handle = net_at(islanded, 0, 0, "handle")
    result = islanded.equilibrium(islanded.capacitance(HANDLE_ALONE), {handle: 300.0})
    (island,) = result.loose
    rows = [k for k, (body, _) in enumerate(result.dofs) if body == island]
    assert not result.displacement[rows].any()
    assert result.force[rows[2]] < 0
    alone = plate.equilibrium(
        plate.capacitance(HANDLE_ALONE), {net_at(plate, 0, 0, "handle"): 300.0}
    )
    assert alone.loose == ()
    assert np.delete(result.displacement, rows) == pytest.approx(alone.displacement, rel=1e-9)
    assert result.softest == pytest.approx(alone.softest)


def test_inputs_that_name_nothing_are_rejected(plate: Structure, shuttle: Structure) -> None:
    model = plate.capacitance(HANDLE_ALONE)
    with pytest.raises(ValueError, match="net 99"):
        plate.equilibrium(model, {99: 1.0})
    with pytest.raises(ValueError, match="another structure"):
        shuttle.equilibrium(model, {})
    with pytest.raises(ValueError, match="max_iterations"):
        StaticOptions(max_iterations=0)
    with pytest.raises(ValueError, match="steps"):
        StaticOptions(steps=0)
    with pytest.raises(ValueError, match="pull_in_resolution"):
        StaticOptions(pull_in_resolution=1.0)

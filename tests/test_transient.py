from dataclasses import replace
from functools import partial

import numpy as np
import pytest
import scipy.integrate
import scipy.interpolate
from numpy.typing import NDArray
from shapely import box

from connflex_core import (
    ElectrostaticParams,
    Layout,
    Material,
    PiecewiseLinear,
    SmallSignalOptions,
    Source,
    Structure,
    TransientOptions,
    extract,
)
from test_equilibrium import (
    BOX,
    EPS0,
    EXTRACTION,
    HANDLE_ALONE,
    SIDE,
    STACK,
    UM,
    net_at,
    plate_on_beams,
    shuttle_layout,
)

PLATES = ElectrostaticParams(dense_cutoff=0.0, plates_cutoff=0.0, interaction_distance=10 * UM)


@pytest.fixture(scope="module")
def plate() -> Structure:
    return plate_on_beams()


@pytest.fixture(scope="module")
def resistive() -> Structure:
    device = replace(STACK.device, material=Material.silicon(resistivity=10.0))
    return extract(shuttle_layout(), replace(STACK, device=device), EXTRACTION)


def plate_mode(plate: Structure) -> tuple[float, float, float]:
    """Spring, mass, and angular frequency of the plate moving toward the handle."""
    modal = plate.modal()
    spring, mass = modal.stiffness[2, 2], modal.mass[2, 2]
    return spring, mass, float(np.sqrt(spring / mass))


def test_source_charges_plate_through_its_resistance(plate: Structure) -> None:
    """One volt through a resistor charges the plate net as an RC circuit: the source
    supplies C V^2, half of it stays in the field and half heats the resistor."""
    model = plate.capacitance(HANDLE_ALONE)
    over = net_at(plate, 50, 50)
    capacitance = model.net_matrix[over, over]
    resistance, volts = 1e6, 1.0
    tau = resistance * capacitance
    pad = plate.electrode_at((-300 * UM, 0.0), "device")
    source = Source(pad, PiecewiseLinear([0.0, 1e-6 * tau], [0.0, volts]), resistance=resistance)
    times = np.linspace(0.0, 5 * tau, 51)
    run = plate.transient(model, [source], times, TransientOptions(quality_factor=10.0))
    assert run.contact is None
    assert run.times == pytest.approx(times)
    rising = volts * (1 - np.exp(-times / tau))
    on_net = model.net == over
    assert run.potentials[:, on_net] == pytest.approx(
        np.repeat(rising[:, None], on_net.sum(), axis=1), abs=1e-5 * volts
    )
    assert run.current[0, 0] == 0.0
    assert run.current[1:, 0] == pytest.approx(
        volts / resistance * np.exp(-times[1:] / tau), abs=1e-5 * volts / resistance
    )
    energy = run.energy
    full = capacitance * volts**2
    assert energy.supplied == pytest.approx(full * (1 - np.exp(-times / tau)), abs=1e-5 * full)
    assert energy.resistive == pytest.approx(
        0.5 * full * (1 - np.exp(-2 * times / tau)), abs=1e-5 * full
    )
    assert np.max(abs(energy.imbalance)) < 1e-6 * full


def test_step_rings_down_to_static_equilibrium(plate: Structure) -> None:
    """A fast step of 30 V sets the plate ringing toward the handle as one damped degree
    of freedom under the pull across its gap, and it settles at the static equilibrium."""
    spring, mass, natural = plate_mode(plate)
    model = plate.capacitance(HANDLE_ALONE)
    over = net_at(plate, 50, 50)
    pad = plate.electrode_at((-300 * UM, 0.0), "device")
    volts, quality, period = 30.0, 5.0, 2 * np.pi / natural
    step = PiecewiseLinear([0.0, 1e-3 * period], [0.0, volts])
    times = np.linspace(0.0, 30 * period, 601)
    options = TransientOptions(quality_factor=quality, tolerance=1e-8)
    run = plate.transient(model, [Source(pad, step)], times, options)
    area = (SIDE * UM) ** 2

    def motion(time: float, state: NDArray[np.float64]) -> list[float]:
        rise, speed = state
        pull = -0.5 * EPS0 * area * step(time) ** 2 / (BOX + rise) ** 2
        return [speed, (pull - spring * rise - mass * natural / quality * speed) / mass]

    exact = scipy.integrate.solve_ivp(
        motion, (0.0, times[-1]), [0.0, 0.0], method="DOP853", t_eval=times, rtol=1e-10, atol=1e-22
    ).y[0]
    settled = -0.5 * EPS0 * area * volts**2 / BOX**2 / spring
    assert run.displacement[:, 2] == pytest.approx(exact, abs=1e-4 * abs(settled))
    rest = plate.equilibrium(model, {over: volts})
    assert run.displacement[-1] == pytest.approx(
        rest.displacement, rel=1e-6, abs=1e-6 * abs(settled)
    )
    assert np.max(abs(run.energy.imbalance)) < 1e-6 * run.energy.supplied[-1]
    assert run.energy.damped[-1] > 0


def test_small_drive_on_bias_follows_small_signal_model(plate: Structure) -> None:
    """A swing of one volt on a bias of 100 V moves the plate as the small-signal
    model predicts once the start has died away."""
    model = plate.capacitance(HANDLE_ALONE)
    over = net_at(plate, 50, 50)
    pad = plate.electrode_at((-300 * UM, 0.0), "device")
    bias, swing, quality = 100.0, 1.0, 5.0
    rest = plate.equilibrium(model, {over: bias})
    signal = plate.small_signal(rest, [pad], SmallSignalOptions(quality_factor=quality))
    mode = int(np.argmax(abs(signal.transduction[:, 0])))
    omega = 0.6 * 2 * np.pi * signal.frequencies_hz[mode]

    def waveform(time: float) -> float:
        return bias + swing * float(np.sin(omega * time))

    period = 2 * np.pi / omega
    times = np.linspace(0.0, 25 * period, 25 * 32 + 1)
    run = plate.transient(
        model,
        [Source(pad, waveform)],
        times,
        TransientOptions(quality_factor=quality, tolerance=1e-8),
        start=rest,
    )
    late = times >= 22 * period
    basis = np.stack(
        [np.cos(omega * times[late]), np.sin(omega * times[late]), np.ones(late.sum())], axis=1
    )
    fit, *_ = np.linalg.lstsq(basis, run.displacement[late, 2], rcond=None)
    measured = fit[0] - 1j * fit[1]
    predicted = signal.displacement([omega / (2 * np.pi)], [-1j * swing])[0, 2]
    assert measured == pytest.approx(predicted, rel=2e-3)
    assert np.max(abs(run.energy.imbalance)) < 1e-6 * np.max(abs(run.energy.supplied))


def test_floating_and_resistive_nets_settle_to_static_equilibrium(resistive: Structure) -> None:
    """Spine driven through its beams, right stator through a resistor, left stator
    floating: from rest the shuttle rings down to the static equilibrium."""
    model = resistive.capacitance(PLATES)
    spine, right, left = (net_at(resistive, x, 0) for x in (0, 250, -250))
    pad = resistive.electrode_at((370 * UM, -240 * UM), "device")
    stator = resistive.electrode_at((250 * UM, 0.0), "device")
    charge = 4e-12
    sources = [Source(pad, 40.0), Source(stator, 10.0, resistance=1e5)]
    rest = resistive.equilibrium(model, {spine: 40.0, right: 10.0}, charges={left: charge})
    slowest = 2 * np.pi * resistive.modal().frequencies_hz[0]
    times = np.linspace(0.0, 60 / slowest * 2 * 3.0, 201)
    run = resistive.transient(
        model,
        sources,
        times,
        TransientOptions(quality_factor=3.0, tolerance=1e-8),
        charges={left: charge},
    )
    assert run.floating == (left,)
    assert run.displacement[-1] == pytest.approx(
        rest.displacement, rel=1e-5, abs=1e-5 * np.max(abs(rest.displacement))
    )
    assert run.potentials[-1] == pytest.approx(rest.potentials, rel=1e-6, abs=1e-6)
    on_left = model.net == left
    assert run.charge[:, on_left].sum(axis=1) == pytest.approx(charge, rel=1e-9)
    assert np.max(abs(run.energy.imbalance)) < 1e-6 * np.max(abs(run.energy.supplied))
    assert run.energy.resistive[-1] > 0


def test_drive_beyond_pull_in_ends_in_contact(plate: Structure) -> None:
    spring, _, natural = plate_mode(plate)
    model = plate.capacitance(HANDLE_ALONE)
    pad = plate.electrode_at((-300 * UM, 0.0), "device")
    pull_in = np.sqrt(8 * spring * BOX**3 / (27 * EPS0 * (SIDE * UM) ** 2))
    period = 2 * np.pi / natural
    step = PiecewiseLinear([0.0, 1e-3 * period], [0.0, 1.3 * pull_in])
    times = np.linspace(0.0, 20 * period, 201)
    run = plate.transient(model, [Source(pad, step)], times, TransientOptions(quality_factor=2.0))
    over = plate.electrode_at((50 * UM, 50 * UM), "device").index
    under = plate.electrode_at((50 * UM, 50 * UM), "handle").index
    assert run.contact == (over, under)
    assert run.times[-1] < times[-1]
    assert run.clearance[-1] == pytest.approx(0.1, rel=1e-6)
    assert np.all(run.clearance[:-1] > 0.1)


PITCH = 30.0
STATORS = {0: (-255, -105, range(-7, -4)), 1: (-65, 85, range(-1, 2)), 2: (125, 275, range(5, 8))}


def toothed_rack() -> Structure:
    """A bar with 21 teeth at a 30 um pitch on two beams 1 mm long, soft along the bar,
    under three anchored stators of three teeth 4 um across, each a third of a pitch on
    from the last."""
    device = [
        box(-330, 0, 330, 40),
        *(box(PITCH * k - 5, 40, PITCH * k + 5, 50) for k in range(-10, 11)),
    ]
    for phase, (low, high, teeth) in STATORS.items():
        device.append(box(low, 64, high, 264))
        device += [
            box(PITCH * m + 10 * phase - 5, 54, PITCH * m + 10 * phase + 5, 64) for m in teeth
        ]
    device += [box(x, -1000, x + 2, 0) for x in (-329, 327)]
    device += [box(-485, -1310, -185, -1000), box(185, -1310, 485, -1000)]
    layout = Layout.from_shapely(device=device, handle=[box(-615, -1400, 615, 400)], unit=UM)
    return extract(layout, STACK, EXTRACTION)


def test_traveling_drive_pulls_out_where_pendulum_reduction_does() -> None:
    """Each stator at V sqrt((1 + cos(w t - 2 pi p / 3)) / 2) makes a field that travels
    one pitch per period and pulls the bar along at a lag. Once its beams and damping ask
    more than the field holds, the bar slips back by a pitch and catches the next tooth.
    The transient on the field method, through the surrogate, follows one degree of
    freedom along the bar under the force of exact field solutions over one pitch. The
    fringe is left out: its rounding where tooth ends pass is finer than the table."""
    rack = toothed_rack()
    modal = rack.modal()
    mass, spring, quality = modal.mass[0, 0], modal.stiffness[0, 0], 0.5
    damping = np.sqrt(spring * mass) / quality
    params = ElectrostaticParams(
        element_size=4 * UM,
        corner_grading=4.0,
        interaction_distance=10 * UM,
        include_plates=False,
        dense_cutoff=float("inf"),
        plates_cutoff=float("inf"),
        surrogate_tolerance=3e-3,
        include_fringe=False,
    )
    model = rack.capacitance(params)
    stators = [
        rack.electrode_at(((low + high) / 2 * UM, 150 * UM), "device")
        for low, high, _ in STATORS.values()
    ]
    along = model.dofs.index((rack.electrode_at((0.0, 20 * UM), "device").body, "ux"))
    volts, omega = 190.0, 2 * np.pi * 400.0

    def level(phase: int, time: float) -> float:
        return float(volts * np.sqrt((1 + np.cos(omega * time - 2 * np.pi * phase / 3)) / 2))

    def levels(time: float) -> NDArray[np.float64]:
        return np.array([level(phase, time) for phase in range(3)])

    picks = np.zeros((3, len(model.net)))
    picks[np.arange(3), [stator.index for stator in stators]] = 1.0
    grid = np.arange(0.0, PITCH, 2.0) * UM
    forces = []
    for travel in grid:
        place = np.zeros(len(model.dofs))
        place[along] = travel
        forces.append(picks @ model.at(place, exact=True).gradient[:, :, along] @ picks.T / 2)
    table = scipy.interpolate.CubicSpline(
        np.append(grid, PITCH * UM), [*forces, forces[0]], axis=0, bc_type="periodic"
    )

    def motion(time: float, state: NDArray[np.float64]) -> list[float]:
        travel, speed = state
        pull = levels(time) @ table(travel) @ levels(time)
        return [speed, (pull - spring * travel - damping * speed) / mass]

    times = np.linspace(0.0, 3e-3, 301)
    reduced = scipy.integrate.solve_ivp(
        motion, (0.0, times[-1]), [0.0, 0.0], method="DOP853", t_eval=times, rtol=1e-9, atol=1e-15
    ).y[0]
    sources = [Source(stator, partial(level, phase)) for phase, stator in enumerate(stators)]
    run = rack.transient(
        model, sources, times, TransientOptions(quality_factor=quality, tolerance=1e-5)
    )
    full = run.displacement[:, run.dofs.index(model.dofs[along])]
    slip = int(np.argmax(reduced))
    assert 0.2 * PITCH * UM < reduced[slip] < 0.8 * PITCH * UM
    assert np.min(reduced[slip:]) < reduced[slip] - PITCH * UM / 2
    assert abs(int(np.argmax(full)) - slip) <= 1
    assert full[:slip] == pytest.approx(reduced[:slip], abs=2e-3 * PITCH * UM)
    assert full == pytest.approx(reduced, abs=2e-2 * PITCH * UM)
    assert np.max(abs(run.energy.imbalance)) < 1e-3 * np.max(run.energy.electrostatic)


def test_waveforms_and_inputs_are_checked(plate: Structure, resistive: Structure) -> None:
    ramp = PiecewiseLinear([0.0, 1.0, 3.0], [0.0, 2.0, 2.0])
    assert ramp(0.5) == 1.0
    assert ramp(-1.0) == 0.0
    assert ramp(5.0) == 2.0
    assert ramp.derivative(0.0) == 2.0
    assert ramp.derivative(1.0) == 0.0
    assert ramp.breakpoints == (0.0, 1.0, 3.0)
    with pytest.raises(ValueError, match="ascending"):
        PiecewiseLinear([0.0, 0.0], [1.0, 2.0])
    with pytest.raises(ValueError, match="resistance"):
        Source(plate.electrodes[0], 1.0, resistance=-1.0)
    with pytest.raises(ValueError, match="contact_clearance"):
        TransientOptions(quality_factor=1.0, contact_clearance=1.0)

    model = plate.capacitance(HANDLE_ALONE)
    options = TransientOptions(quality_factor=10.0)
    pad = plate.electrode_at((-300 * UM, 0.0), "device")
    middle = plate.electrode_at((50 * UM, 50 * UM), "device")
    with pytest.raises(TypeError, match="Source"):
        plate.transient(model, [pad], [0.0, 1e-6], options)  # type: ignore[list-item]
    with pytest.raises(TypeError, match="numbered apart"):
        plate.transient(model, [Source(3, 1.0)], [0.0, 1e-6], options)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="one node"):
        plate.transient(model, [Source(pad, 1.0), Source(middle, 1.0)], [0.0, 1e-6], options)
    with pytest.raises(ValueError, match="floating"):
        plate.transient(model, [Source(pad, 1.0)], [0.0, 1e-6], options, charges={pad.net: 0.0})
    with pytest.raises(ValueError, match="times"):
        plate.transient(model, [Source(pad, 1.0)], [1e-6, 0.0], options)
    with pytest.raises(ValueError, match="another structure"):
        resistive.transient(model, [], [0.0, 1e-6], options)

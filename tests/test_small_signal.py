import shutil
import subprocess
import warnings
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
import scipy.linalg
import scipy.optimize
from numpy.typing import NDArray

from connflex_core import (
    Electrode,
    ElectrostaticParams,
    Equilibrium,
    LinearizationWarning,
    Material,
    SmallSignal,
    SmallSignalOptions,
    Structure,
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

VOLTS = 300.0
QUALITY = 1000.0
PLATES = ElectrostaticParams(dense_cutoff=0.0, plates_cutoff=0.0, interaction_distance=10 * UM)


@pytest.fixture(scope="module")
def plate() -> Structure:
    return plate_on_beams()


@pytest.fixture(scope="module")
def resistive() -> Structure:
    """Comb shuttle whose beams have resistance, so its spine lies behind them."""
    device = replace(STACK.device, material=Material.silicon(resistivity=10.0))
    return extract(shuttle_layout(), replace(STACK, device=device), EXTRACTION)


def plate_closed_form(plate: Structure) -> dict[str, float]:
    """Spring, mass, gap, softened spring, and force per volt of the plate at VOLTS."""
    modal = plate.modal()
    spring, mass, area = modal.stiffness[2, 2], modal.mass[2, 2], (SIDE * UM) ** 2

    def unbalance(rise: float) -> float:
        return float(spring * rise + 0.5 * EPS0 * area * VOLTS**2 / (BOX + rise) ** 2)

    gap = BOX + scipy.optimize.brentq(unbalance, -BOX / 3, 0.0, xtol=1e-18)
    return {
        "mass": mass,
        "gap": gap,
        "softened": spring - EPS0 * area * VOLTS**2 / gap**3,
        "slope": -EPS0 * area * VOLTS / gap**2,
    }


def plate_signal(plate: Structure) -> tuple[Equilibrium, SmallSignal]:
    rest = plate.equilibrium(plate.capacitance(HANDLE_ALONE), {net_at(plate, 50, 50): VOLTS})
    pad = plate.electrode_at((-300 * UM, 0.0), "device")
    return rest, plate.small_signal(rest, [pad], SmallSignalOptions(quality_factor=QUALITY))


def test_plate_resonator_matches_closed_form_motional_branch(plate: Structure) -> None:
    """Plate over the handle at a bias: one mode moves it, softened by the field, and
    seen from its pad it is a capacitor in parallel with a series RLC branch."""
    exact = plate_closed_form(plate)
    rest, signal = plate_signal(plate)
    mode = int(np.argmax(abs(signal.transduction[:, 0])))
    natural = np.sqrt(exact["softened"] / exact["mass"])
    assert signal.frequencies_hz[mode] == pytest.approx(natural / (2 * np.pi), rel=1e-8)
    assert abs(signal.modes[2, mode]) == pytest.approx(1 / np.sqrt(exact["mass"]), rel=1e-8)
    unbiased = np.sqrt(plate.modal().stiffness[2, 2] / exact["mass"]) / (2 * np.pi)
    assert signal.frequencies_hz[mode] < unbiased
    assert abs(signal.transduction[mode, 0]) == pytest.approx(
        abs(exact["slope"]) / np.sqrt(exact["mass"]), rel=1e-8
    )

    inductance, capacitance, resistance = signal.motional()
    assert inductance[mode, 0] == pytest.approx(exact["mass"] / exact["slope"] ** 2, rel=1e-8)
    assert capacitance[mode, 0] == pytest.approx(exact["slope"] ** 2 / exact["softened"], rel=1e-8)
    assert resistance[mode, 0] == pytest.approx(
        np.sqrt(exact["softened"] * exact["mass"]) / (QUALITY * exact["slope"] ** 2), rel=1e-8
    )

    over = net_at(plate, 50, 50)
    static = rest.model.net_matrix[over, over]
    assert signal.capacitance == pytest.approx(np.array([[static]]), rel=1e-12)
    hertz = signal.frequencies_hz[mode] * np.array([0.5, 0.999, 1.0, 1.001, 2.0])
    omega = 2 * np.pi * hertz
    branch = exact["slope"] ** 2 / exact["mass"]
    expected = (
        1j * omega * (static + branch / (natural**2 - omega**2 + 1j * omega * natural / QUALITY))
    )
    assert signal.admittance(hertz)[:, 0, 0] == pytest.approx(expected, rel=1e-6)


def biased_shuttle(resistive: Structure) -> tuple[Equilibrium, list[Electrode]]:
    """Spine at 40 V through its beams, left stator floating, right stator and one
    spine pad as ports, handle at zero."""
    model = resistive.capacitance(PLATES)
    spine, left = net_at(resistive, 0, 0), net_at(resistive, -250, 0)
    rest = resistive.equilibrium(model, {spine: 40.0}, charges={left: 4e-12})
    stator = resistive.electrode_at((250 * UM, 0.0), "device")
    pad = resistive.electrode_at((370 * UM, -240 * UM), "device")
    return rest, [stator, pad]


def direct_admittance(
    structure: Structure,
    rest: Equilibrium,
    electrodes: list[Electrode],
    hertz: NDArray[np.float64],
) -> NDArray[np.complex128]:
    """Port admittance from one linear solve per frequency over every degree of
    freedom and every electrode, floating ones included, with modal damping of Q."""
    modal = structure.modal()
    mass, model, volts = modal.mass, rest.model, rest.potentials
    ports = [electrode.index for electrode in electrodes]
    count = len(structure.electrodes)
    pull = np.einsum("abk,b->ak", model.gradient, volts)
    spring = modal.stiffness + model.stiffness(volts)
    squared, shapes = scipy.linalg.eigh(rest.stiffness, mass)
    damping = mass @ shapes @ np.diag(np.sqrt(squared) / QUALITY) @ shapes.T @ mass
    conductance = np.zeros((count, count))
    for link in structure.links:
        a, b = link.electrodes
        assert link.resistance
        conductance[np.ix_([a, b], [a, b])] += np.array([[1, -1], [-1, 1]]) / link.resistance
    driven = {model.net[p] for p in ports} | set(rest.floating)
    held = list(ports) + [a for a in range(count) if model.net[a] not in driven]
    dofs = len(rest.dofs)
    result = np.zeros((len(hertz), len(ports), len(ports)), dtype=complex)
    for f, omega in enumerate(2 * np.pi * hertz):
        system = np.block(
            [
                [spring - omega**2 * mass + 1j * omega * damping, -pull.T],
                [1j * omega * pull, conductance + 1j * omega * model.matrix],
            ]
        )
        for column, port in enumerate(ports):
            rows = system.copy()
            given = np.zeros(dofs + count, dtype=complex)
            for a in held:
                rows[dofs + a] = 0.0
                rows[dofs + a, dofs + a] = 1.0
            given[dofs + port] = 1.0
            state = np.linalg.solve(rows, given)
            current = system[dofs:] @ state
            result[f, :, column] = current[ports]
    return result


def test_admittance_matches_direct_solve_of_network(resistive: Structure) -> None:
    rest, ports = biased_shuttle(resistive)
    signal = resistive.small_signal(rest, ports, SmallSignalOptions(quality_factor=QUALITY))
    assert signal.ports == tuple(ports)
    spine = {e.index for e in resistive.electrodes if e.net == net_at(resistive, 0, 0)}
    assert signal.nodes[:2] == ((ports[0].index,), (ports[1].index,))
    assert {node for (node,) in signal.nodes[2:]} == spine - {ports[1].index}
    mode = int(np.argmax(abs(signal.transduction[:, 0])))
    hertz = np.concatenate(
        [[1e3, 3e4], signal.frequencies_hz[mode] * np.array([0.5, 0.999, 1.0, 1.001, 2.0])]
    )
    admittance = signal.admittance(hertz)
    assert admittance == pytest.approx(direct_admittance(resistive, rest, ports, hertz), rel=1e-6)
    assert admittance == pytest.approx(admittance.transpose(0, 2, 1), rel=1e-9)


def element_values(text: str) -> dict[str, list[str]]:
    return {line.split()[0]: line.split() for line in text.splitlines() if line[:1].isalpha()}


def test_netlist_of_plate_is_capacitor_beside_series_branch(plate: Structure) -> None:
    exact = plate_closed_form(plate)
    _, signal = plate_signal(plate)
    mode = int(np.argmax(abs(signal.transduction[:, 0])))
    text = signal.netlist("plate", modes=[mode])
    lines = [line for line in text.splitlines() if line and not line.startswith("*")]
    assert lines[0] == ".subckt plate port0"
    assert lines[-1] == ".ends plate"
    elements = element_values(text)
    assert set(elements) == {f"{kind}{mode}" for kind in ("Lm", "Rm", "Cm", "Vm")} | {
        "Cg0",
        f"Em{mode}_0",
        f"Fm{mode}_0",
    }
    assert float(elements["Cg0"][3]) == pytest.approx(signal.capacitance[0, 0], rel=1e-10)
    assert float(elements[f"Lm{mode}"][3]) == pytest.approx(
        exact["mass"] / exact["slope"] ** 2, rel=1e-8
    )
    assert float(elements[f"Cm{mode}"][3]) == pytest.approx(
        exact["slope"] ** 2 / exact["softened"], rel=1e-8
    )
    assert abs(float(elements[f"Em{mode}_0"][5])) == pytest.approx(1.0)
    assert float(elements[f"Fm{mode}_0"][4]) == float(elements[f"Em{mode}_0"][5])


@pytest.mark.skipif(shutil.which("ngspice") is None, reason="ngspice is not installed")
def test_netlist_reproduces_admittance_in_ngspice(resistive: Structure, tmp_path: Path) -> None:
    rest, ports = biased_shuttle(resistive)
    signal = resistive.small_signal(rest, ports, SmallSignalOptions(quality_factor=QUALITY))
    low, high = 0.5 * signal.frequencies_hz[0], 2.0 * signal.frequencies_hz[0]
    hertz = np.linspace(low, high, 7)
    computed = signal.admittance(hertz)
    for column in range(len(ports)):
        sources = "\n".join(
            f"V{p} port{p} 0 DC 0 AC {1 if p == column else 0}" for p in range(len(ports))
        )
        deck = (
            f"* admittance column {column}\n{signal.netlist('device')}\n"
            f"X1 {' '.join(f'port{p}' for p in range(len(ports)))} device\n{sources}\n"
            f".control\nset wr_singlescale\noption numdgt=15\n"
            f"ac lin 7 {float(low):.15g} {float(high):.15g}\n"
            f"wrdata out.txt {' '.join(f'i(V{p})' for p in range(len(ports)))}\n"
            f".endc\n.end\n"
        )
        (tmp_path / "deck.cir").write_text(deck)
        subprocess.run(["ngspice", "-b", "deck.cir"], cwd=tmp_path, check=True, capture_output=True)
        table = np.loadtxt(tmp_path / "out.txt")
        assert table[:, 0] == pytest.approx(hertz, rel=1e-9)
        drawn = -(table[:, 1::2] + 1j * table[:, 2::2])
        assert drawn == pytest.approx(computed[:, :, column], rel=1e-6)


def test_linearity_compares_drive_with_bias_and_motion_with_gap(plate: Structure) -> None:
    exact = plate_closed_form(plate)
    _, signal = plate_signal(plate)
    mode = int(np.argmax(abs(signal.transduction[:, 0])))
    natural = np.sqrt(exact["softened"] / exact["mass"])
    hertz = 0.5 * signal.frequencies_hz[mode]
    omega = 2 * np.pi * hertz
    swing = 0.01
    travel = abs(
        exact["slope"]
        * swing
        / (exact["mass"] * (natural**2 - omega**2 + 1j * omega * natural / QUALITY))
    )
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        small = signal.linearity([swing], hertz)
    assert small.voltage == pytest.approx(swing / VOLTS, rel=1e-9)
    assert small.gap == pytest.approx(travel / exact["gap"], rel=1e-6)
    assert abs(small.peak[2]) == pytest.approx(travel, rel=1e-8)
    with pytest.warns(LinearizationWarning, match="gap"):
        large = signal.linearity([20.0], signal.frequencies_hz[mode])
    assert large.gap > 0.1


def test_inputs_outside_small_signal_are_rejected(plate: Structure, resistive: Structure) -> None:
    rest, signal = plate_signal(plate)
    options = SmallSignalOptions(quality_factor=QUALITY)
    pad = plate.electrode_at((-300 * UM, 0.0), "device")
    middle = plate.electrode_at((50 * UM, 50 * UM), "device")
    with pytest.raises(ValueError, match="one node"):
        plate.small_signal(rest, [pad, middle], options)
    with pytest.raises(TypeError, match="numbered apart"):
        plate.small_signal(rest, [net_at(plate, 50, 50)], options)  # type: ignore[list-item]
    with pytest.raises(ValueError, match="electrode 0 is of another structure"):
        plate.small_signal(rest, [resistive.electrodes[0]], options)
    with pytest.raises(ValueError, match="another structure"):
        resistive.small_signal(rest, [pad], options)
    with pytest.raises(ValueError, match="quality_factor"):
        SmallSignalOptions(quality_factor=0.0)
    with pytest.raises(ValueError, match="frequencies"):
        signal.admittance([-1.0])
    with pytest.raises(ValueError, match="drive"):
        signal.linearity([1.0, 2.0], 1e3)

    shuttle, _ = biased_shuttle(resistive)
    left = resistive.electrode_at((-250 * UM, 0.0), "device")
    with pytest.raises(ValueError, match="floating"):
        resistive.small_signal(shuttle, [left], options)

    area, spring = (SIDE * UM) ** 2, plate.modal().stiffness[2, 2]

    def unbalance(rise: float) -> float:
        return float(spring * rise + 0.5 * EPS0 * area * VOLTS**2 / (BOX + rise) ** 2)

    tipping = scipy.optimize.brentq(unbalance, -0.999 * BOX, -BOX / 3, xtol=1e-18)
    start = np.zeros(len(rest.dofs))
    start[2] = 1.02 * tipping
    model = plate.capacitance(HANDLE_ALONE)
    unstable = plate.equilibrium(model, {net_at(plate, 50, 50): VOLTS}, start=start)
    assert not unstable.stable
    with pytest.raises(ValueError, match="stable"):
        plate.small_signal(unstable, [pad], options)

from dataclasses import replace
from itertools import pairwise

import numpy as np
import pytest
from shapely import box

from connflex_core import (
    CapacitanceModel,
    ElectrostaticParams,
    Layout,
    Structure,
    extract,
)
from connflex_core.electrostatics.dense import charge_derivatives, solve_charges
from connflex_core.electrostatics.outlines import Mesh
from connflex_core.progress import quiet
from test_equilibrium import EXTRACTION, STACK, net_at

UM = 1e-6
PITCH = 16.0


def elements(points: list[tuple[float, float]], size: float) -> list[tuple[complex, complex]]:
    pieces = []
    for a, b in pairwise(points):
        start, end = complex(*a), complex(*b)
        count = max(1, int(np.ceil(abs(end - start) / size)))
        fractions = np.linspace(0.0, 1.0, count + 1)
        pieces += [
            (start + (end - start) * fractions[i], start + (end - start) * fractions[i + 1])
            for i in range(count)
        ]
    return pieces


def comb_cell(x: float, size: float) -> tuple[list[tuple[complex, complex]], ...]:
    """Stator bar on top with a finger down, shuttle bar below with a finger up, 4 um
    fingers and gaps, one pitch wide starting at ``x``, in micrometers."""
    stator = elements([(x + 16, 50), (x, 50)], size) + elements(
        [(x, 40), (x + 2, 40), (x + 2, 0), (x + 6, 0), (x + 6, 40), (x + 16, 40)], size
    )
    shuttle = elements([(x, -20), (x + 16, -20)], size) + elements(
        [(x + 16, -10), (x + 14, -10), (x + 14, 30), (x + 10, 30), (x + 10, -10), (x, -10)],
        size,
    )
    return stator, shuttle


def mesh_of(conductors: list[list[tuple[complex, complex]]], shift: complex = 0) -> Mesh:
    start = np.array([a for pieces in conductors for a, _ in pieces]) + shift
    end = np.array([b for pieces in conductors for _, b in pieces]) + shift
    owner = np.array([k for k, pieces in enumerate(conductors) for _ in pieces], dtype=np.intp)
    return Mesh(
        np.column_stack([start.real, start.imag]) * UM,
        np.column_stack([end.real, end.imag]) * UM,
        owner,
    )


def test_periodic_cell_is_the_limit_of_long_finite_arrays() -> None:
    """The middle cell of a finite comb approaches the cell of the endless comb as the
    comb grows, as one over its length."""
    period = np.array([PITCH * UM, 0.0])
    cell = solve_charges(mesh_of(list(comb_cell(0.0, 2.0))), 2, quiet, period=period)
    endless = -cell.capacitance[0, 1]
    assert cell.capacitance[0, 1] == pytest.approx(cell.capacitance[1, 0], rel=1e-10)
    middle = {}
    for count in (9, 17):
        stator: list[tuple[complex, complex]] = []
        shuttle: list[tuple[complex, complex]] = []
        for k in range(count):
            top, bottom = comb_cell(k * PITCH, 2.0)
            stator += top
            shuttle += bottom
        width = count * PITCH
        stator += elements([(width, 40), (width, 50)], 2.0) + elements([(0, 50), (0, 40)], 2.0)
        shuttle += elements([(0, -10), (0, -20)], 2.0) + elements([(width, -20), (width, -10)], 2.0)
        mesh = mesh_of([stator, shuttle])
        center = (mesh.start[:, 0] + mesh.end[:, 0]) / 2 / UM
        low = (count // 2) * PITCH
        core = (center >= low) & (center < low + PITCH)
        middle[count] = -solve_charges(mesh, 2, quiet, core).capacitance[0, 1]
    assert abs(middle[17] - endless) < abs(middle[9] - endless)
    limit = (17 * middle[17] - 9 * middle[9]) / 8
    assert limit == pytest.approx(endless, rel=1e-3)


def test_periodic_derivatives_match_finite_differences() -> None:
    """Travel of the shuttle in every cell at once keeps the row periodic, and the
    derivatives of the periodic kernel follow it exactly."""
    period = np.array([PITCH * UM, 0.0])
    stator, shuttle = comb_cell(0.0, 2.0)
    center = complex(8.0, 15.0)
    pivots = np.array([[center.real, center.imag]] * 2) * UM

    def capacitance(shift: complex = 0) -> float:
        moved = [(a + shift, b + shift) for a, b in shuttle]
        solved = solve_charges(mesh_of([stator, moved]), 2, quiet, period=period)
        return float(solved.capacitance[0, 1])

    charges = solve_charges(mesh_of([stator, shuttle]), 2, quiet, period=period)
    gradient, hessian = charge_derivatives(charges, np.array([0, 1]), pivots, quiet)
    step = 0.02
    along = (capacitance(complex(0, step)) - capacitance(complex(0, -step))) / (2 * step * UM)
    assert gradient[0, 1, 1, 1] == pytest.approx(along, rel=1e-4)
    across = (capacitance(complex(step, 0)) - capacitance(complex(-step, 0))) / (2 * step * UM)
    assert gradient[0, 1, 1, 0] == pytest.approx(across, abs=1e-4 * abs(along))
    curve = (capacitance(complex(step, 0)) - 2 * capacitance() + capacitance(complex(-step, 0))) / (
        step * UM
    ) ** 2
    assert hessian[0, 1, 1, 0, 1, 0] == pytest.approx(curve, rel=1e-3)
    assert gradient[0, 1, 0, :2] == pytest.approx(-gradient[0, 1, 1, :2], abs=1e-6 * abs(along))


WINDOWS = ElectrostaticParams(
    element_size=4 * UM,
    corner_grading=4.0,
    interaction_distance=10 * UM,
    include_plates=False,
    dense_cutoff=0.0,
    plates_cutoff=float("inf"),
    window_elements=200,
)


@pytest.fixture(scope="module")
def long_comb() -> Structure:
    """Forty fingers on a released bar between two beams, against forty fingers of an
    anchored stator, 4 um fingers and gaps at a 16 um pitch."""
    device = [box(-340, 0, 340, 40), box(-340, 94, 340, 254)]
    device += [box(16 * k + 10, 40, 16 * k + 14, 80) for k in range(-20, 20)]
    device += [box(16 * k + 2, 54, 16 * k + 6, 94) for k in range(-20, 20)]
    device += [box(x, -300, x + 10, 0) for x in (-340, 330)]
    device += [box(-490, -610, -190, -300), box(190, -610, 490, -300)]
    layout = Layout.from_shapely(device=device, handle=[box(-600, -700, 600, 400)], unit=UM)
    return extract(layout, STACK, EXTRACTION)


def test_long_comb_is_found_as_one_row(long_comb: Structure) -> None:
    model = long_comb.capacitance(replace(WINDOWS, periodic_min_cells=8))
    rows = [row for problem in model.rest.groups[0].plan.problems for row in problem.rows]
    assert len(rows) == 1
    assert abs(rows[0].period[0]) == pytest.approx(16 * UM)
    assert rows[0].period[1] == pytest.approx(0.0, abs=1e-12)
    assert 30 <= len(rows[0].cells) <= 38


def test_periodic_row_matches_whole_solution_at_rest_and_displaced(long_comb: Structure) -> None:
    """The row takes the capacitance of its interior from one cell, and the derivatives
    in turn from a window around it. Travel along and across the fingers and a turn of
    the shuttle agree with one dense solution of the whole comb, at rest and with the
    shuttle moved and turned."""
    whole = long_comb.capacitance(replace(WINDOWS, dense_cutoff=float("inf")))
    periodic = long_comb.capacitance(replace(WINDOWS, periodic_min_cells=8))
    stator = net_at(long_comb, 0, 200)
    shuttle = net_at(long_comb, 0, 20)
    turn = periodic.dofs.index((1, "tz"))
    potentials = np.where(periodic.net == stator, 50.0, 0.0)
    for place in (np.zeros(len(periodic.dofs)), _moved(periodic, turn)):
        truth, found = whole.at(place), periodic.at(place)
        assert found.net_matrix[stator, shuttle] == pytest.approx(
            truth.net_matrix[stator, shuttle], rel=4e-3
        )
        pick = np.ix_(model_net(found, stator), model_net(found, shuttle))

        def slope(model: CapacitanceModel, k: int, pick: tuple[np.ndarray, ...] = pick) -> float:
            return float(model.gradient[pick][..., k].sum())

        along = slope(truth, 1)
        assert slope(found, 1) == pytest.approx(along, rel=2e-3)
        assert slope(found, 0) == pytest.approx(slope(truth, 0), abs=2e-3 * abs(along))
        assert slope(found, turn) == pytest.approx(
            slope(truth, turn), abs=2e-3 * abs(along) * 300 * UM
        )
        for k in (0, turn):
            assert found.stiffness(potentials)[k, k] == pytest.approx(
                truth.stiffness(potentials)[k, k], rel=5e-3
            )


def _moved(model: CapacitanceModel, turn: int) -> np.ndarray:
    place = np.zeros(len(model.dofs))
    place[1] = 0.5 * UM
    place[0] = 0.3 * UM
    place[turn] = 1e-4
    return place


def model_net(model: CapacitanceModel, net: int) -> np.ndarray:
    return np.flatnonzero(model.net == net)

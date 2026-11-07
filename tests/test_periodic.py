from dataclasses import replace
from itertools import pairwise
from typing import Any, cast

import numpy as np
import pytest
from shapely import Polygon, affinity, box, union_all

from connflex_core import (
    CapacitanceModel,
    ElectrostaticParams,
    ExtractionParams,
    Layout,
    Structure,
    extract,
    sidewall_capacitance,
)
from connflex_core.electrostatics import periodic, windows
from connflex_core.electrostatics.dense import Charges, charge_derivatives, solve_charges
from connflex_core.electrostatics.outlines import Mesh
from connflex_core.electrostatics.periodic import _remainder
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
    model = long_comb.capacitance(WINDOWS)
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
    periodic = long_comb.capacitance(WINDOWS)
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


def two_sided_comb(offset: tuple[float, float]) -> Structure:
    """Twenty fingers on either side of a released spine, each side against its own
    anchored stator, 1.5 um fingers and gaps at a 6 um pitch, drawn ``offset`` from the
    origin in micrometers; the tips on the two sides are mirror images across the spine."""
    pitch, finger, length, overlap, spine = 6.0, 1.5, 40.0, 30.0, 12.0
    count, half = 20, 60.0
    device = [box(-half - 12, -spine / 2, half + 12, spine / 2)]
    for sign in (1, -1):
        root = sign * spine / 2
        tip = root + sign * (length - overlap)
        for k in range(count):
            x = -half + pitch * k
            low, high = sorted((root, root + sign * length))
            device.append(box(x, low, x + finger, high))
        for k in range(-1, count):
            x = -half + pitch * k + finger + 1.5
            low, high = sorted((tip, tip + sign * length))
            device.append(box(x, low, x + finger, high))
        base = tip + sign * length
        low, high = sorted((base, base + sign * 60))
        device.append(box(-half - 32, low, half + 32, high))
    for sign in (1, -1):
        x1, x2 = sorted((sign * (half + 12), sign * (half + 162)))
        device.append(box(x1, -2, x2, 2))
        p1, p2 = sorted((sign * (half + 162), sign * (half + 222)))
        device.append(box(p1, -30, p2, 30))
    device = [affinity.translate(part, *offset) for part in device]
    handle = affinity.translate(box(-400, -300, 400, 300), *offset)
    layout = Layout.from_shapely(device=device, handle=[handle], unit=UM)
    undercut = ExtractionParams(
        release_hole_max_area=100 * UM * UM, release_distance=7 * UM, body_min_width=10 * UM
    )
    return extract(layout, STACK, undercut)


def rows_of(model: CapacitanceModel) -> list[int]:
    return [
        len(row.cells)
        for group in model.rest.groups
        for problem in group.plan.problems
        for row in problem.rows
    ]


def test_mirror_image_tips_do_not_break_a_row(monkeypatch: pytest.MonkeyPatch) -> None:
    """Elements of a cell are ordered by position; the tips of fingers on the two sides of
    a spine share a position along the row up to rounding, which must not reorder them
    from one cell to the next, wherever the comb is drawn."""
    params = ElectrostaticParams(
        element_size=0.75 * UM,
        interaction_distance=20 * UM,
        include_plates=False,
        include_fringe=False,
    )
    found = [
        rows_of(two_sided_comb(offset).capacitance(params))
        for offset in ((0.0, 0.0), (3000.0, 3000.0))
    ]
    assert found[0] == found[1]
    assert found[0]
    assert min(found[0]) >= 8


def test_remainder_stays_finite_far_across_the_row(monkeypatch: pytest.MonkeyPatch) -> None:
    """Across the row the periodic kernel grows without bound while sine overflows;
    beyond a few periods the remainder takes its asymptotic form, which agrees with
    the direct one where both can be evaluated."""
    reach = PITCH * UM
    for sign in (1.0, -1.0):
        point = np.array([0.3 + sign * 18j], dtype=np.complex128) * reach / np.pi
        direct = _remainder(point, reach)
        monkeypatch.setattr(periodic, "_FAR", 10.0)
        asymptotic = _remainder(point, reach)
        monkeypatch.setattr(periodic, "_FAR", 20.0)
        for close, beyond in zip(direct, asymptotic, strict=True):
            assert beyond[0] == pytest.approx(close[0], rel=1e-10)
        apart = np.array([0.3 + sign * 800j], dtype=np.complex128) * reach / np.pi
        value, slope, curve = _remainder(apart, reach)
        assert np.isfinite(value).all()
        assert np.isfinite(slope).all()
        assert np.isfinite(curve).all()
        assert value[0] > 700
        assert slope[0] == pytest.approx(-sign * 1j * np.pi / reach, rel=1e-2)
        assert abs(curve[0]) < 1e-3 * (np.pi / reach) ** 2


def test_row_evaluated_far_across_itself_stays_finite(long_comb: Structure) -> None:
    periodic = long_comb.capacitance(WINDOWS)
    stator, shuttle = net_at(long_comb, 0, 200), net_at(long_comb, 0, 20)
    place = np.zeros(len(periodic.dofs))
    place[periodic.dofs.index((1, "uy"))] = -120 * PITCH * UM
    found = periodic.at(place)
    assert np.isfinite(found.matrix).all()
    assert np.isfinite(found.gradient).all()
    assert abs(found.net_matrix[stator, shuttle]) < 1e-2 * abs(periodic.net_matrix[stator, shuttle])


def test_windows_beside_translated_rows_are_solved_once(monkeypatch: pytest.MonkeyPatch) -> None:
    """Three combs 150 um apart: what each row leaves to the windows, bar ends and the
    walls behind its fingers, is cut in the frame of the row, so the windows of every comb
    take the solutions of the first comb's."""

    def comb(y: float) -> list[Polygon]:
        shuttle = [
            box(-160, 0, 160, 15),
            *(box(PITCH * k + 10, 15, PITCH * k + 14, 55) for k in range(-10, 10)),
        ]
        stator = [
            box(-160, 69, 160, 84),
            *(box(PITCH * k + 2, 29, PITCH * k + 6, 69) for k in range(-10, 10)),
        ]
        return [
            cast(
                Polygon,
                affinity.scale(affinity.translate(union_all(part), 0, y), UM, UM, origin=(0, 0)),
            )
            for part in (shuttle, stator)
        ]

    params = ElectrostaticParams(
        element_size=4 * UM,
        corner_grading=4.0,
        interaction_distance=15 * UM,
        dense_cutoff=0.0,
        plates_cutoff=float("inf"),
        window_elements=100,
    )
    solved = []

    def counted(*args: Any, **kwargs: Any) -> Charges:
        solved.append(1)
        return solve_charges(*args, **kwargs)

    monkeypatch.setattr(windows, "solve_charges", counted)
    one = sidewall_capacitance(comb(0), params)
    once = len(solved)
    three = sidewall_capacitance([*comb(0), *comb(150), *comb(300)], params)
    assert len(solved) == 2 * once
    assert three.windows == 3 * one.windows
    for k in range(3):
        mine = np.arange(2 * k, 2 * k + 2)
        assert three.capacitance[np.ix_(mine, mine)] == pytest.approx(one.capacitance, rel=1e-12)
        assert three.gradient[np.ix_(mine, mine, mine)] == pytest.approx(
            one.gradient, abs=1e-12 * np.max(abs(one.gradient))
        )


def _moved(model: CapacitanceModel, turn: int) -> np.ndarray:
    place = np.zeros(len(model.dofs))
    place[1] = 0.5 * UM
    place[0] = 0.3 * UM
    place[turn] = 1e-4
    return place


def model_net(model: CapacitanceModel, net: int) -> np.ndarray:
    return np.flatnonzero(model.net == net)

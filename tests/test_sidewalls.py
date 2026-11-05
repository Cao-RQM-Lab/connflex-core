import warnings

import numpy as np
import pytest
from shapely import MultiPolygon, Polygon, box

from connflex_core import ElectrostaticParams, sidewall_capacitance
from connflex_core.electrostatics.outlines import mesh_outlines

EPS0 = 8.8541878188e-12
UM = 1e-6


def circle(radius: float, sides: int, x: float = 0.0, y: float = 0.0) -> Polygon:
    angle = 2 * np.pi * (np.arange(sides) + 0.5) / sides
    return Polygon(np.column_stack([x + radius * np.cos(angle), y + radius * np.sin(angle)]))


def um(polygon: Polygon) -> Polygon:
    exterior = np.asarray(polygon.exterior.coords) * UM
    return Polygon(exterior, [np.asarray(ring.coords) * UM for ring in polygon.interiors])


def params(size: float, grading: float = 16.0, limit: int = 8000) -> ElectrostaticParams:
    return ElectrostaticParams(element_size=size * UM, corner_grading=grading, max_elements=limit)


def test_eccentric_circles_have_the_closed_form_capacitance() -> None:
    a, b, offset = 1.0, 2.5, 0.6
    inner = um(circle(a, 180, x=offset))
    outer = um(Polygon(circle(4.0, 180).exterior, [circle(b, 180).exterior]))
    result = sidewall_capacitance([inner, outer], params(1.0))
    expected = 2 * np.pi * EPS0 / np.arccosh((a * a + b * b - offset * offset) / (2 * a * b))
    c = result.capacitance
    assert c[0, 0] == pytest.approx(expected, rel=1e-4)
    assert c == pytest.approx(np.array([[1, -1], [-1, 1]]) * c[0, 0], rel=1e-12)
    outside = np.hypot(*(result.mesh.start / UM).T) > 3.0
    assert outside.sum() == 180
    assert abs(result.charge[outside, 0].sum()) < 1e-5 * c[0, 0]


def test_long_parallel_strips_approach_the_parallel_plate() -> None:
    gap = 2.0

    def strips(length: float) -> list[Polygon]:
        return [um(box(0, 0, length, 4)), um(box(0, 4 + gap, length, 8 + gap))]

    short = sidewall_capacitance(strips(80.0), params(1.0))
    long = sidewall_capacitance(strips(160.0), params(1.0))
    added = long.capacitance[0, 0] - short.capacitance[0, 0]
    assert added == pytest.approx(EPS0 * 80.0 / gap, rel=1e-2)
    assert long.capacitance[0, 0] > EPS0 * 160.0 / gap

    middle = (long.mesh.start + long.mesh.end) / 2 / UM
    facing = np.argmin(np.hypot(middle[:, 0] - 80.0, middle[:, 1] - 4.0))
    length = np.hypot(*(long.mesh.end[facing] - long.mesh.start[facing]))
    assert long.mesh.conductor[facing] == 0
    assert long.charge[facing, 0] / length == pytest.approx(EPS0 / (gap * UM), rel=1e-6)


def test_refining_the_mesh_converges_from_below() -> None:
    squares = [um(box(0, 0, 10, 10)), um(box(12, 0, 22, 10))]
    coarse, medium, fine = (
        sidewall_capacitance(squares, params(size)).capacitance[0, 0] for size in (2.0, 1.0, 0.5)
    )
    assert coarse < medium < fine
    assert fine - medium < 0.5 * (medium - coarse)
    uniform = sidewall_capacitance(squares, params(1.0, grading=1.0)).capacitance[0, 0]
    assert fine - medium < 0.2 * (fine - uniform)


def test_cavity_that_holds_no_conductor_is_not_meshed() -> None:
    solid, other = um(box(0, 0, 10, 10)), um(box(12, 0, 22, 10))
    hollow = Polygon(solid.exterior, [um(box(3, 3, 7, 7)).exterior])
    filled = sidewall_capacitance([solid, other], params(1.0))
    perforated = sidewall_capacitance([hollow, other], params(1.0))
    assert len(perforated.mesh.conductor) == len(filled.mesh.conductor)
    assert np.array_equal(perforated.capacitance, filled.capacitance)


def test_elements_tile_a_long_outline_within_the_size_limit() -> None:
    strip = MultiPolygon([um(box(0, 0, 3000, 7))])
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        mesh = mesh_outlines([strip], 2.0 * UM, 16.0)
    lengths = np.hypot(*(mesh.end - mesh.start).T)
    assert lengths.sum() == pytest.approx(strip.length)
    assert lengths.max() <= 2.0 * UM
    assert lengths.min() == pytest.approx(2.0 * UM / 16 * np.expm1(1.0), rel=0.2)


def test_inputs_that_cannot_be_solved_are_rejected() -> None:
    squares = [um(box(0, 0, 10, 10)), um(box(12, 0, 22, 10))]
    with pytest.raises(ValueError, match="element_size"):
        sidewall_capacitance(squares, ElectrostaticParams())
    with pytest.raises(ValueError, match="max_elements"):
        sidewall_capacitance(squares, params(1.0, limit=50))
    with pytest.raises(ValueError, match="conductors 0 and 1"):
        sidewall_capacitance([squares[0], um(box(10, 0, 20, 10))], params(1.0))

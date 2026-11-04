import numpy as np
import pytest
from shapely import MultiPolygon, Polygon, box, unary_union

from connflex_core import Stack
from connflex_core.geometry.heal import as_multipolygon
from connflex_core.mechanics.inertia import AreaMoments, area_moments, mass_properties

UM = 1e-6


def moments_tuple(moments: AreaMoments) -> list[float]:
    return [moments.area, moments.x, moments.y, moments.xx, moments.yy, moments.xy]


def test_rectangle_moments_are_exact() -> None:
    moments = area_moments(MultiPolygon([box(1, 2, 4, 6)]))
    assert moments_tuple(moments) == pytest.approx([12, 30, 48, 84, 208, 120])


def test_moments_add_over_concave_and_holed_shapes() -> None:
    arm_x, arm_y = box(0, 0, 4, 1), box(0, 1, 1, 4)
    ell = area_moments(as_multipolygon(unary_union([arm_x, arm_y])))
    parts = np.add(
        moments_tuple(area_moments(MultiPolygon([arm_x]))),
        moments_tuple(area_moments(MultiPolygon([arm_y]))),
    )
    assert moments_tuple(ell) == pytest.approx(list(parts))

    outer, inner = box(0, 0, 10, 10), box(4, 4, 6, 6)
    ring = area_moments(MultiPolygon([Polygon(outer.exterior, [inner.exterior])]))
    difference = np.subtract(
        moments_tuple(area_moments(MultiPolygon([outer]))),
        moments_tuple(area_moments(MultiPolygon([inner]))),
    )
    assert moments_tuple(ring) == pytest.approx(list(difference))


def test_moments_ignore_ring_orientation() -> None:
    outer = [(0, 0), (10, 0), (10, 10), (0, 10)]
    inner = [(4, 4), (6, 4), (6, 6), (4, 6)]
    drawn = MultiPolygon([Polygon(outer, [inner])])
    reversed_rings = MultiPolygon([Polygon(outer[::-1], [inner[::-1]])])
    assert moments_tuple(area_moments(drawn)) == pytest.approx([96, 480, 480, 3232, 3232, 2400])
    assert moments_tuple(area_moments(reversed_rings)) == pytest.approx(
        moments_tuple(area_moments(drawn))
    )


def test_single_plate_matches_the_cuboid_formulas(stack: Stack) -> None:
    a, b, t = 100 * UM, 50 * UM, stack.device.thickness
    x0, y0 = 4000 * UM, -3000 * UM
    plate = MultiPolygon([box(x0, y0, x0 + a, y0 + b)])
    props = mass_properties([plate, MultiPolygon(), MultiPolygon()], stack)
    mass = 2330 * t * a * b
    assert props.mass == pytest.approx(mass)
    assert props.center == pytest.approx([x0 + a / 2, y0 + b / 2, -t / 2])
    expected = mass / 12 * np.diag([b * b + t * t, a * a + t * t, a * a + b * b])
    assert props.inertia == pytest.approx(expected, abs=1e-11 * np.max(expected))
    assert props.matrix[:3, :3] == pytest.approx(mass * np.eye(3))
    assert props.matrix[3:, 3:] == pytest.approx(props.inertia)
    assert props.matrix[:3, 3:] == pytest.approx(np.zeros((3, 3)))


def test_two_layers_shift_the_center_and_add_parallel_axis_terms(stack: Stack) -> None:
    a, b = 100 * UM, 50 * UM
    t_d, t_b = stack.device.thickness, stack.box.thickness
    plate = MultiPolygon([box(0, 0, a, b)])
    props = mass_properties([plate, plate, MultiPolygon()], stack)
    m_d, m_b = 2330 * t_d * a * b, 2200 * t_b * a * b
    z_d, z_b = -t_d / 2, -t_d - t_b / 2
    cz = (m_d * z_d + m_b * z_b) / (m_d + m_b)
    assert props.center[2] == pytest.approx(cz)
    ixx = sum(
        m * ((b * b + t * t) / 12 + (z - cz) ** 2) for m, t, z in [(m_d, t_d, z_d), (m_b, t_b, z_b)]
    )
    assert props.inertia[0, 0] == pytest.approx(ixx)
    assert props.inertia[0, 2] == pytest.approx(0.0, abs=1e-9 * ixx)

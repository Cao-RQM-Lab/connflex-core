import numpy as np
import pytest
from shapely import box

from connflex_core import Beam, Material
from connflex_core.mechanics.stiffness import (
    effective_stiffness,
    effective_torsion_constant,
    local_stiffness,
    torsion_constant,
)

UM = 1e-6
SILICON = Material.silicon()


def beam(length: float, width: float, thickness: float, angle: float = 0.0) -> Beam:
    end = (length * np.cos(angle), length * np.sin(angle), -thickness / 2)
    return Beam(
        index=0,
        bodies=(0, 1),
        length=length,
        width=width,
        thickness=thickness,
        start=(0.0, 0.0, -thickness / 2),
        end=end,
        polygon=box(0, -width / 2, length, width / 2),
        net=0,
    )


def test_torsion_constant_limits() -> None:
    assert torsion_constant(100.0, 1.0) == pytest.approx(100 / 3, rel=1e-2)
    assert torsion_constant(1.0, 1.0) == pytest.approx(0.1406, rel=2e-3)
    assert torsion_constant(1.0, 3.0) == torsion_constant(3.0, 1.0)


def test_short_beam_correction_vanishes_for_zero_beta_and_long_beams() -> None:
    j = torsion_constant(2.0, 1.0)
    assert effective_torsion_constant(2.0, 1.0, 10.0, beta=0.0) == j
    assert effective_torsion_constant(2.0, 1.0, 2000.0, beta=0.63) == pytest.approx(j, rel=1e-5)
    assert effective_torsion_constant(2.0, 1.0, 2.0, beta=0.63) == pytest.approx(1.63 * j)


def test_local_stiffness_is_symmetric_with_six_rigid_body_modes() -> None:
    length = 200 * UM
    k = local_stiffness(beam(length, 10 * UM, 25 * UM), SILICON)
    assert k == pytest.approx(k.T)
    scale = np.diag([1, 1, 1, 1 / length, 1 / length, 1 / length] * 2)
    assert np.linalg.matrix_rank(scale @ k @ scale) == 6
    assert np.all(np.linalg.eigvalsh(scale @ k @ scale) > -1e-9 * np.max(np.abs(k)))
    modes = []
    for axis in range(3):
        translation = np.zeros(12)
        translation[axis] = translation[6 + axis] = 1.0
        rotation = np.zeros(12)
        rotation[3 + axis] = rotation[9 + axis] = 1.0
        rotation[6:9] = np.cross(np.eye(3)[axis], [length, 0.0, 0.0])
        modes += [translation, rotation]
    for mode in modes:
        assert np.linalg.norm(k @ mode) <= 1e-9 * np.linalg.norm(k) * np.linalg.norm(mode)


def test_cantilever_deflections_match_closed_forms() -> None:
    length, width, thickness = 200 * UM, 10 * UM, 25 * UM
    b = beam(length, width, thickness)
    anchor_center, tip_center = np.array([-50 * UM, 0.0, -200 * UM]), np.array(b.end)
    k_tip = effective_stiffness(b, SILICON, anchor_center, tip_center)[6:, 6:]
    e = SILICON.youngs_modulus
    for axis, stiffness in [
        (0, e * width * thickness / length),
        (1, e * thickness * width**3 / (4 * length**3)),
        (2, e * width * thickness**3 / (4 * length**3)),
    ]:
        force = np.zeros(6)
        force[axis] = 1e-6
        displacement = np.linalg.solve(k_tip, force)
        assert displacement[axis] == pytest.approx(force[axis] / stiffness, rel=1e-9)


def test_effective_stiffness_is_invariant_under_rotation() -> None:
    straight = beam(200 * UM, 10 * UM, 25 * UM)
    rotated = beam(200 * UM, 10 * UM, 25 * UM, angle=0.7)
    offset_a, offset_b = np.array([-30 * UM, 20 * UM, -100 * UM]), np.array(
        [40 * UM, -10 * UM, 0.0]
    )
    c, s = np.cos(0.7), np.sin(0.7)
    rotate = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    k_straight = effective_stiffness(straight, SILICON, offset_a, np.array(straight.end) + offset_b)
    k_rotated = effective_stiffness(
        rotated, SILICON, rotate @ offset_a, np.array(rotated.end) + rotate @ offset_b
    )
    assert k_straight == pytest.approx(k_straight.T)
    reference = np.linalg.eigvalsh(k_straight)
    assert np.linalg.eigvalsh(k_rotated) == pytest.approx(
        reference, rel=1e-9, abs=1e-12 * np.max(reference)
    )


def test_lever_arm_maps_the_beam_end_stiffness_to_an_offset_center() -> None:
    """A body whose center lies off the beam end sees the end stiffness through the lever
    of equations (27) and (28): a motion of the center carries the end by the cross
    product of the turn with the arm."""
    b = beam(200 * UM, 10 * UM, 25 * UM)
    anchor_center = np.array([-50 * UM, 0.0, -200 * UM])
    arm = np.array([30 * UM, -20 * UM, -40 * UM])
    at_end = effective_stiffness(b, SILICON, anchor_center, np.array(b.end))
    offset = effective_stiffness(b, SILICON, anchor_center, np.array(b.end) + arm)
    skew = np.array([[0.0, -arm[2], arm[1]], [arm[2], 0.0, -arm[0]], [-arm[1], arm[0], 0.0]])
    carry = np.eye(12)
    carry[6:9, 9:12] = skew  # end travel = center travel + turn x (end - center)
    expected = carry.T @ at_end @ carry
    assert offset == pytest.approx(expected, rel=1e-12, abs=1e-12 * np.max(abs(at_end)))
    assert offset[6:, 6:] != pytest.approx(at_end[6:, 6:])

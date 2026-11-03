"""Beam stiffness on body degrees of freedom."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from connflex_core.materials import Material
from connflex_core.model import Beam

Matrix = NDArray[np.float64]


def torsion_constant(width: float, thickness: float) -> float:
    """Saint-Venant torsion constant of a solid rectangle, equation (20)."""
    a, b = max(width, thickness), min(width, thickness)
    return a * b**3 * (1 / 3 - 0.21 * (b / a) * (1 - b**4 / (12 * a**4)))


def effective_torsion_constant(width: float, thickness: float, length: float, beta: float) -> float:
    """Torsion constant with short-beam stiffening factor, equation (21)."""
    a = max(width, thickness)
    return torsion_constant(width, thickness) * (1 + beta * a**2 / length**2)


def local_stiffness(beam: Beam, material: Material) -> Matrix:
    """12 by 12 Euler-Bernoulli stiffness in the beam frame, equations (16) to (19).

    Local x runs along the beam, y lies in plane, and z points out of plane. Degrees of
    freedom are ordered as in equation (12) for each end.
    """
    e, g = material.youngs_modulus, material.shear_modulus
    if g is None:
        raise ValueError("material shear modulus is not resolved")
    length, w, t = beam.length, beam.width, beam.thickness
    j = effective_torsion_constant(w, t, length, material.torsion_beta)
    k_a, k_t = e * w * t / length, g * j / length
    k_y, c_z, g_z = (
        e * t * w**3 / length**3,
        e * t * w**3 / (2 * length**2),
        e * t * w**3 / (3 * length),
    )
    k_z, c_y, g_y = (
        e * w * t**3 / length**3,
        e * w * t**3 / (2 * length**2),
        e * w * t**3 / (3 * length),
    )
    k00 = np.array(
        [
            [k_a, 0, 0, 0, 0, 0],
            [0, k_y, 0, 0, 0, c_z],
            [0, 0, k_z, 0, -c_y, 0],
            [0, 0, 0, k_t, 0, 0],
            [0, 0, -c_y, 0, g_y, 0],
            [0, c_z, 0, 0, 0, g_z],
        ]
    )
    k11 = np.array(
        [
            [k_a, 0, 0, 0, 0, 0],
            [0, k_y, 0, 0, 0, -c_z],
            [0, 0, k_z, 0, c_y, 0],
            [0, 0, 0, k_t, 0, 0],
            [0, 0, c_y, 0, g_y, 0],
            [0, -c_z, 0, 0, 0, g_z],
        ]
    )
    k01 = np.array(
        [
            [-k_a, 0, 0, 0, 0, 0],
            [0, -k_y, 0, 0, 0, c_z],
            [0, 0, -k_z, 0, -c_y, 0],
            [0, 0, 0, -k_t, 0, 0],
            [0, 0, c_y, 0, g_y / 2, 0],
            [0, -c_z, 0, 0, 0, g_z / 2],
        ]
    )
    return np.block([[k00, k01], [k01.T, k11]])


def transformation(beam: Beam, center_a: Matrix, center_b: Matrix) -> Matrix:
    """12 by 12 map from beam-end forces to body forces, equation (31).

    Each block is the end offset from the body center of mass followed by rotation from
    beam frame to global frame.
    """
    start, end = np.asarray(beam.start), np.asarray(beam.end)
    direction = end - start
    alpha = float(np.arctan2(direction[1], direction[0]))
    rotation = _rotation6(alpha)
    t_a = _offset(start - center_a) @ rotation
    t_b = _offset(end - center_b) @ rotation
    return np.block([[t_a, np.zeros((6, 6))], [np.zeros((6, 6)), t_b]])


def effective_stiffness(
    beam: Beam, material: Material, center_a: Matrix, center_b: Matrix
) -> Matrix:
    """12 by 12 stiffness felt by both bodies, equation (33)."""
    t = transformation(beam, center_a, center_b)
    return t @ local_stiffness(beam, material) @ t.T


def _rotation6(alpha: float) -> Matrix:
    c, s = np.cos(alpha), np.sin(alpha)
    r3 = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    r6 = np.zeros((6, 6))
    r6[:3, :3] = r3
    r6[3:, 3:] = r3
    return r6


def _offset(d: Matrix) -> Matrix:
    skew = np.array([[0.0, -d[2], d[1]], [d[2], 0.0, -d[0]], [-d[1], d[0], 0.0]])
    s = np.eye(6)
    s[3:, :3] = skew
    return s

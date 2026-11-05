"""Integrals of the logarithmic kernel over pairs of straight elements."""

from __future__ import annotations

from math import comb

import numpy as np
from numpy.typing import NDArray

Array = NDArray[np.float64]
Complex = NDArray[np.complex128]

SERIES_ORDER = 7
SERIES_REACH = 0.125
_EVEN = [
    [comb(2 * m, 2 * k) / ((2 * k + 1) * (2 * m - 2 * k + 1)) for k in range(m + 1)]
    for m in range(SERIES_ORDER + 1)
]


def pair_integrals(a0: Complex, a1: Complex, b0: Complex, b1: Complex) -> Array:
    """Double integral of ln|r - r'| over element pairs, equation (64).

    Points are complex numbers x + iy; ``a0`` and ``a1`` are ends of the first element of
    each pair, ``b0`` and ``b1`` of the second, broadcast against each other. Pairs
    farther apart than their lengths by ``1 / SERIES_REACH`` use the series (69), the
    others the closed form (67). A pair may share an end point or be one element twice,
    but two different elements must not cross or overlap.
    """
    a0, a1, b0, b1 = np.broadcast_arrays(a0, a1, b0, b1)
    center = (a0 + a1 - b0 - b1) / 2
    near = np.abs(a1 - a0) + np.abs(b1 - b0) > 2 * SERIES_REACH * np.abs(center)
    integral = _series(a0, a1, b0, b1, np.where(near, 1.0, center))
    index = np.nonzero(near)
    integral[index] = _closed_form(a0[index], a1[index], b0[index], b1[index])
    return integral


def _closed_form(a0: Complex, a1: Complex, b0: Complex, b1: Complex) -> Array:
    """Second difference of G over the four pairs of end points, equations (66) and (67)."""
    tangents = (a1 - a0) / np.abs(a1 - a0) * (b1 - b0) / np.abs(b1 - b0)
    center = (a0 + a1 - b0 - b1) / 2
    same = center == 0
    turn = np.where(same, 1.0, np.conj(center) / np.where(same, 1.0, np.abs(center)))
    total = np.zeros(a0.shape, dtype=np.complex128)
    for u, sign in ((a1 - b1, 1.0), (a1 - b0, -1.0), (a0 - b1, -1.0), (a0 - b0, 1.0)):
        apart = u != 0
        log = np.log(np.where(apart, u * turn, 1.0))
        total += sign * np.where(apart, u * u * (log / 2 - 0.75), 0.0)
    integral: Array = (-total / tangents).real
    return integral


def _series(a0: Complex, a1: Complex, b0: Complex, b1: Complex, center: Complex) -> Array:
    """Expansion about the difference of element centers, equation (69)."""
    first, second = a1 - a0, b1 - b0
    alpha, beta = (first / (2 * center)) ** 2, (second / (2 * center)) ** 2
    alphas, betas = [np.ones_like(alpha)], [np.ones_like(beta)]
    for _ in range(SERIES_ORDER):
        alphas.append(alphas[-1] * alpha)
        betas.append(betas[-1] * beta)
    series = np.zeros(alpha.shape, dtype=np.complex128)
    for m in range(1, SERIES_ORDER + 1):
        for k in range(m + 1):
            series += (_EVEN[m][k] / (2 * m)) * alphas[k] * betas[m - k]
    integral: Array = np.abs(first) * np.abs(second) * (np.log(np.abs(center)) - series.real)
    return integral

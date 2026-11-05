"""Integrals of the logarithmic kernel over pairs of straight elements."""

from __future__ import annotations

from math import comb
from typing import overload

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
_ODD = [
    [0.0]
    + [comb(2 * m - 1, 2 * k - 1) / ((2 * k + 1) * (2 * m - 2 * k + 1)) for k in range(1, m + 1)]
    for m in range(SERIES_ORDER + 1)
]


@overload
def pair_integrals(
    a0: Complex, a1: Complex, b0: Complex, b1: Complex, pivot: None = None
) -> tuple[Array, None]: ...


@overload
def pair_integrals(
    a0: Complex, a1: Complex, b0: Complex, b1: Complex, pivot: Complex
) -> tuple[Array, Array]: ...


def pair_integrals(
    a0: Complex, a1: Complex, b0: Complex, b1: Complex, pivot: Complex | None = None
) -> tuple[Array, Array | None]:
    """Double integral of ln|r - r'| over element pairs, equation (64), and its slopes.

    Points are complex numbers x + iy; ``a0`` and ``a1`` are ends of the first element of
    each pair, ``b0`` and ``b1`` of the second, broadcast against each other. Pairs
    farther apart than their lengths by ``1 / SERIES_REACH`` use the series (69), the
    others the closed form (67). A pair may share an end point or be one element twice,
    but two different elements must not cross or overlap.

    With ``pivot``, the second result stacks the derivatives of the integral for travel
    of the first element along x, along y, and about ``pivot``, equation (78).
    """
    a0, a1, b0, b1 = np.broadcast_arrays(a0, a1, b0, b1)
    center = (a0 + a1 - b0 - b1) / 2
    near = np.abs(a1 - a0) + np.abs(b1 - b0) > 2 * SERIES_REACH * np.abs(center)
    turning = None if pivot is None else np.broadcast_to(pivot, a0.shape)
    integral, slopes = _series(a0, a1, b0, b1, np.where(near, 1.0, center), turning)
    index = np.nonzero(near)
    close = _closed_form(
        a0[index], a1[index], b0[index], b1[index], None if turning is None else turning[index]
    )
    integral[index] = close[0]
    if slopes is not None and close[1] is not None:
        for slope, fixed in zip(slopes, close[1], strict=True):
            slope[index] = fixed
    return integral, slopes


def _closed_form(
    a0: Complex, a1: Complex, b0: Complex, b1: Complex, pivot: Complex | None
) -> tuple[Array, Array | None]:
    """Second differences over the four pairs of end points, equations (67) and (79)."""
    tangents = (a1 - a0) / np.abs(a1 - a0) * (b1 - b0) / np.abs(b1 - b0)
    center = (a0 + a1 - b0 - b1) / 2
    same = center == 0
    turn = np.where(same, 1.0, np.conj(center) / np.where(same, 1.0, np.abs(center)))
    total = np.zeros(a0.shape, dtype=np.complex128)
    force = np.zeros(a0.shape, dtype=np.complex128)
    moment = np.zeros(a0.shape, dtype=np.complex128)
    for u, sign, end in (
        (a1 - b1, 1.0, a1),
        (a1 - b0, -1.0, a1),
        (a0 - b1, -1.0, a0),
        (a0 - b0, 1.0, a0),
    ):
        apart = u != 0
        log = np.log(np.where(apart, u * turn, 1.0))
        total += sign * np.where(apart, u * u * (log / 2 - 0.75), 0.0)
        if pivot is not None:
            first = sign * np.where(apart, u * (log - 1), 0.0)
            force += first
            moment += first * (end - pivot)
    integral = -total / tangents
    if pivot is None:
        return integral.real, None
    force, moment = -force / tangents, -moment / tangents - integral
    return integral.real, np.stack([force.real, -force.imag, -moment.imag])


def _series(
    a0: Complex, a1: Complex, b0: Complex, b1: Complex, center: Complex, pivot: Complex | None
) -> tuple[Array, Array | None]:
    """Expansion about the difference of element centers, equations (69) and (80)."""
    first, second = a1 - a0, b1 - b0
    alpha, beta = (first / (2 * center)) ** 2, (second / (2 * center)) ** 2
    alphas, betas = [np.ones_like(alpha)], [np.ones_like(beta)]
    for _ in range(SERIES_ORDER):
        alphas.append(alphas[-1] * alpha)
        betas.append(betas[-1] * beta)
    series = np.zeros(alpha.shape, dtype=np.complex128)
    even = np.zeros(alpha.shape, dtype=np.complex128)
    odd = np.zeros(alpha.shape, dtype=np.complex128)
    for m in range(1, SERIES_ORDER + 1):
        for k in range(m + 1):
            term = alphas[k] * betas[m - k]
            series += (_EVEN[m][k] / (2 * m)) * term
            if pivot is not None:
                even += _EVEN[m][k] * term
                odd += _ODD[m][k] * term
    area = np.abs(first) * np.abs(second)
    integral: Array = area * (np.log(np.abs(center)) - series.real)
    if pivot is None:
        return integral, None
    force = area / center * (1 + even)
    moment = ((a0 + a1) / 2 - pivot) * force - area * odd
    return integral, np.stack([force.real, -force.imag, -moment.imag])

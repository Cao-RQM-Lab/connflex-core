"""Integrals of the logarithmic kernel over pairs of straight elements."""

from __future__ import annotations

from collections.abc import Iterator
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
_TURNS = np.array([0.0, 0.0, 1.0, 0.0, 0.0, 1.0])
_BASIS = np.zeros((6, 3), dtype=np.complex128)
_BASIS[[0, 1, 2], [0, 0, 1]] = 1, 1j, 1j
_BASIS[[3, 4, 5], [0, 0, 2]] = -1, -1j, -1j


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
    alpha, beta, safe = _ratios(a0, a1, b0, b1, np.where(near, 1.0, center))
    series = np.zeros(a0.shape, dtype=np.complex128)
    for m, _, term in _terms(alpha, beta):
        series += term / (2 * m)
    integral: Array = np.abs(a1 - a0) * np.abs(b1 - b0) * (np.log(np.abs(safe)) - series.real)
    integral[near] = _corner_sums(a0[near], a1[near], b0[near], b1[near])[0].real
    return integral


def pair_derivatives(
    a0: Complex, a1: Complex, b0: Complex, b1: Complex, pivot_a: Complex, pivot_b: Complex
) -> tuple[Array, Array]:
    """Gradient and Hessian of the pair integral for rigid motion of both elements.

    The six parameters are travel along x and y and rotation about ``pivot_a`` of the
    first element, then the same about ``pivot_b`` of the second, equations (83) to (88).
    The first axis of the gradient and the first two of the Hessian run over them. The
    two elements of a pair must not touch.
    """
    a0, a1, b0, b1, pivot_a, pivot_b = np.broadcast_arrays(a0, a1, b0, b1, pivot_a, pivot_b)
    center = (a0 + a1 - b0 - b1) / 2
    near = np.abs(a1 - a0) + np.abs(b1 - b0) > 2 * SERIES_REACH * np.abs(center)
    gradient, hessian = _series_derivatives(
        a0, a1, b0, b1, np.where(near, 1.0, center), pivot_a, pivot_b
    )
    if near.any():
        ends = (a0[near], a1[near], b0[near], b1[near])
        close = _closed_derivatives(*ends, pivot_a[near], pivot_b[near])
        gradient[:, near] = close[0]
        hessian[:, :, near] = close[1]
    return gradient, hessian


def _corner_sums(
    a0: Complex,
    a1: Complex,
    b0: Complex,
    b1: Complex,
    weights: tuple[Complex, Complex] | None = None,
) -> tuple[Complex, Complex, Complex]:
    """Second differences of G, of G' times weights, and of G'' times pairs of weights.

    Weights are 1, the first end point less ``weights[0]``, and the second end point less
    ``weights[1]``. The sums carry the factor -1 / (tau_i tau_j) of equation (67).
    """
    scale = -1 / ((a1 - a0) / np.abs(a1 - a0) * (b1 - b0) / np.abs(b1 - b0))
    center = (a0 + a1 - b0 - b1) / 2
    same = center == 0
    turn = np.where(same, 1.0, np.conj(center) / np.where(same, 1.0, np.abs(center)))
    value = np.zeros(a0.shape, dtype=np.complex128)
    first = np.zeros((3, *a0.shape), dtype=np.complex128)
    second = np.zeros((3, 3, *a0.shape), dtype=np.complex128)
    for u, sign, end_a, end_b in (
        (a1 - b1, 1.0, a1, b1),
        (a1 - b0, -1.0, a1, b0),
        (a0 - b1, -1.0, a0, b1),
        (a0 - b0, 1.0, a0, b0),
    ):
        apart = u != 0
        log = np.log(np.where(apart, u * turn, 1.0))
        value += sign * np.where(apart, u * u * (log / 2 - 0.75), 0.0)
        if weights is None:
            continue
        basis = (np.ones(a0.shape, dtype=np.complex128), end_a - weights[0], end_b - weights[1])
        slope = sign * u * (log - 1)
        for i in range(3):
            first[i] += slope * basis[i]
            for j in range(i, 3):
                second[i, j] += sign * log * basis[i] * basis[j]
                second[j, i] = second[i, j]
    return scale * value, scale * first, scale * second


def _closed_derivatives(
    a0: Complex, a1: Complex, b0: Complex, b1: Complex, pivot_a: Complex, pivot_b: Complex
) -> tuple[Array, Array]:
    """Derivatives of the closed form (67), equations (85) and (86)."""
    value, first, second = _corner_sums(a0, a1, b0, b1, (pivot_a, pivot_b))
    turns = _TURNS.reshape(6, *([1] * a0.ndim))
    pull = np.einsum("pa,a...->p...", _BASIS, first)
    gradient = pull - 1j * turns * value
    hessian = np.einsum("pa,qb,ab...->pq...", _BASIS, _BASIS, second)
    hessian -= 1j * (turns[:, None] * pull[None, :] + pull[:, None] * turns[None, :])
    hessian -= turns[:, None] * turns[None, :] * value
    hessian[2, 2] -= first[1]
    hessian[5, 5] += first[2]
    return gradient.real, hessian.real


def _series_derivatives(
    a0: Complex,
    a1: Complex,
    b0: Complex,
    b1: Complex,
    center: Complex,
    pivot_a: Complex,
    pivot_b: Complex,
) -> tuple[Array, Array]:
    """Derivatives of the series (69), equations (87) and (88)."""
    alpha, beta, center = _ratios(a0, a1, b0, b1, center)
    sums = np.zeros((9, *a0.shape), dtype=np.complex128)
    for m, k, term in _terms(alpha, beta):
        rest = m - k
        for row, weight in enumerate(
            (
                1,
                2 * m + 1,
                k / m,
                rest / m,
                k,
                rest,
                2 * k * k / m,
                2 * k * rest / m,
                2 * rest * rest / m,
            )
        ):
            if weight:
                sums[row] += weight * term
    slope = (1 + sums[0]) / center
    bend = -(1 + sums[1]) / (center * center)
    arm_a, arm_b = (a0 + a1) / 2 - pivot_a, (b0 + b1) / 2 - pivot_b
    one = np.ones(a0.shape, dtype=np.complex128)
    travel = np.stack([one, 1j * one, 1j * arm_a, -one, -1j * one, -1j * arm_b])
    gradient = slope * travel
    gradient[2] -= 1j * sums[2]
    gradient[5] -= 1j * sums[3]
    twist = np.zeros((6, *a0.shape), dtype=np.complex128)
    twist[2], twist[5] = 2j * sums[4] / center, 2j * sums[5] / center
    hessian = bend * travel[:, None] * travel[None, :]
    hessian += twist[:, None] * travel[None, :] + travel[:, None] * twist[None, :]
    hessian[2, 2] += sums[6] - slope * arm_a
    hessian[5, 5] += sums[8] + slope * arm_b
    hessian[2, 5] += sums[7]
    hessian[5, 2] += sums[7]
    area = np.abs(a1 - a0) * np.abs(b1 - b0)
    return area * gradient.real, area * hessian.real


def _ratios(
    a0: Complex, a1: Complex, b0: Complex, b1: Complex, center: Complex
) -> tuple[Complex, Complex, Complex]:
    """Squares of half element vectors over the difference of centers, as in (69)."""
    return ((a1 - a0) / (2 * center)) ** 2, ((b1 - b0) / (2 * center)) ** 2, center


def _terms(alpha: Complex, beta: Complex) -> Iterator[tuple[int, int, Complex]]:
    """Terms of S_m in (69), each with its m and k."""
    alphas, betas = [np.ones_like(alpha)], [np.ones_like(beta)]
    for _ in range(SERIES_ORDER):
        alphas.append(alphas[-1] * alpha)
        betas.append(betas[-1] * beta)
    for m in range(1, SERIES_ORDER + 1):
        for k in range(m + 1):
            yield m, k, _EVEN[m][k] * alphas[k] * betas[m - k]

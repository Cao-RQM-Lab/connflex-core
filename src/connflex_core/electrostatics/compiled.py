"""Pair integrals and derivatives of ``kernel`` compiled with numba, one pair at a time.

Each function takes the arguments of its numpy namesake, broadcasts them, and runs the
pairs in one compiled loop that releases the interpreter lock, so that threads share it.
"""

from __future__ import annotations

import cmath
import math

import numba
import numpy as np
from numpy.typing import NDArray

from connflex_core.electrostatics.kernel import _EVEN, SERIES_ORDER, SERIES_REACH

Array = NDArray[np.float64]
Complex = NDArray[np.complex128]

ORDER = SERIES_ORDER
REACH = SERIES_REACH
EVEN = np.zeros((ORDER + 1, ORDER + 1))
for _m, _row in enumerate(_EVEN):
    EVEN[_m, : len(_row)] = _row
HALVES = np.array([0.0] + [0.5 / m for m in range(1, ORDER + 1)])
# Weights of the eight sums beyond the first of (87) and (88) carried by each term.
WEIGHTS = np.zeros((ORDER + 1, ORDER + 1, 8))
for _m in range(1, ORDER + 1):
    for _k in range(_m + 1):
        _rest = _m - _k
        WEIGHTS[_m, _k] = (
            2 * _m + 1,
            _k / _m,
            _rest / _m,
            _k,
            _rest,
            2 * _k * _k / _m,
            2 * _k * _rest / _m,
            2 * _rest * _rest / _m,
        )
# Corner sum each parameter of (85) pulls on, with its factor, and which parameters turn.
SLOT = np.array([0, 0, 1, 0, 0, 2])
FACTOR = np.array([1.0, 1j, 1j, -1.0, -1j, -1j])
TURNS = np.array([0.0, 0.0, 1.0, 0.0, 0.0, 1.0])

# No NaN, infinity, or signed zero, and fused multiply-add: without the first two LLVM calls
# a library routine for every complex product.
FLAGS = {"nnan", "ninf", "nsz", "contract"}
_inner = numba.njit(cache=True, nogil=True, inline="always", fastmath=FLAGS)


@_inner
def _powers(alpha, beta, ap, bp):
    ap[0] = 1.0
    bp[0] = 1.0
    for k in range(ORDER):
        ap[k + 1] = ap[k] * alpha
        bp[k + 1] = bp[k] * beta


@_inner
def _corner(u, turn):
    """Term of a corner in the closed form (67): zero where the corner is shared."""
    if u == 0:
        return 0j
    return u * u * (cmath.log(u * turn) / 2 - 0.75)


@_inner
def _integral(a0, a1, b0, b1, center, ap, bp):
    la, lb = abs(a1 - a0), abs(b1 - b0)
    if la + lb > 2 * REACH * abs(center):
        scale = -1 / ((a1 - a0) / la * (b1 - b0) / lb)
        turn = 1.0 + 0j if center == 0 else center.conjugate() / abs(center)
        value = _corner(a1 - b1, turn) - _corner(a1 - b0, turn)
        value += _corner(a0 - b0, turn) - _corner(a0 - b1, turn)
        return (scale * value).real
    d = (a1 - a0) / (2 * center)
    e = (b1 - b0) / (2 * center)
    _powers(d * d, e * e, ap, bp)
    series = 0j
    for m in range(1, ORDER + 1):
        term = 0j
        for k in range(m + 1):
            term += EVEN[m, k] * ap[k] * bp[m - k]
        series += term * HALVES[m]
    return la * lb * (math.log(abs(center)) - series.real)


@numba.njit(cache=True, nogil=True, fastmath=FLAGS)
def _integrals(a0, a1, b0, b1, out):
    ap = np.empty(ORDER + 1, np.complex128)
    bp = np.empty(ORDER + 1, np.complex128)
    for n in range(out.size):
        center = (a0[n] + a1[n] - b0[n] - b1[n]) / 2
        out[n] = _integral(a0[n], a1[n], b0[n], b1[n], center, ap, bp)


@_inner
def _closed(a0, a1, b0, b1, center, pa, pb, first, second, pull, gradient, hessian, n):
    """Column ``n`` of the derivatives by the closed form, equations (85) and (86)."""
    scale = -1 / ((a1 - a0) / abs(a1 - a0) * (b1 - b0) / abs(b1 - b0))
    turn = 1.0 + 0j if center == 0 else center.conjugate() / abs(center)
    value = 0j
    for i in range(3):
        first[i] = 0j
        for j in range(3):
            second[i, j] = 0j
    for corner in range(4):
        if corner == 0:
            u, sign, ea, eb = a1 - b1, 1.0, a1, b1
        elif corner == 1:
            u, sign, ea, eb = a1 - b0, -1.0, a1, b0
        elif corner == 2:
            u, sign, ea, eb = a0 - b1, -1.0, a0, b1
        else:
            u, sign, ea, eb = a0 - b0, 1.0, a0, b0
        if u == 0:
            continue
        log = cmath.log(u * turn)
        value += sign * u * u * (log / 2 - 0.75)
        slope = sign * u * (log - 1)
        weight = sign * log
        arm_a, arm_b = ea - pa, eb - pb
        first[0] += slope
        first[1] += slope * arm_a
        first[2] += slope * arm_b
        second[0, 0] += weight
        second[0, 1] += weight * arm_a
        second[0, 2] += weight * arm_b
        second[1, 1] += weight * arm_a * arm_a
        second[1, 2] += weight * arm_a * arm_b
        second[2, 2] += weight * arm_b * arm_b
    value *= scale
    for i in range(3):
        first[i] *= scale
        for j in range(i, 3):
            second[i, j] *= scale
            second[j, i] = second[i, j]
    for p in range(6):
        pull[p] = FACTOR[p] * first[SLOT[p]]
        gradient[p, n] = (pull[p] - 1j * TURNS[p] * value).real
    for p in range(6):
        for q in range(6):
            entry = FACTOR[p] * FACTOR[q] * second[SLOT[p], SLOT[q]]
            entry -= 1j * (TURNS[p] * pull[q] + pull[p] * TURNS[q])
            entry -= TURNS[p] * TURNS[q] * value
            hessian[p, q, n] = entry.real
    hessian[2, 2, n] -= first[1].real
    hessian[5, 5, n] += first[2].real


@_inner
def _series(a0, a1, b0, b1, center, pa, pb, ap, bp, travel, twist, gradient, hessian, n):
    """Column ``n`` of the derivatives by the series, equations (87) and (88)."""
    d = (a1 - a0) / (2 * center)
    e = (b1 - b0) / (2 * center)
    _powers(d * d, e * e, ap, bp)
    s0 = s1 = s2 = s3 = s4 = s5 = s6 = s7 = s8 = 0j
    for m in range(1, ORDER + 1):
        for k in range(m + 1):
            term = EVEN[m, k] * ap[k] * bp[m - k]
            s0 += term
            s1 += WEIGHTS[m, k, 0] * term
            s2 += WEIGHTS[m, k, 1] * term
            s3 += WEIGHTS[m, k, 2] * term
            s4 += WEIGHTS[m, k, 3] * term
            s5 += WEIGHTS[m, k, 4] * term
            s6 += WEIGHTS[m, k, 5] * term
            s7 += WEIGHTS[m, k, 6] * term
            s8 += WEIGHTS[m, k, 7] * term
    slope = (1 + s0) / center
    bend = -(1 + s1) / (center * center)
    arm_a = (a0 + a1) / 2 - pa
    arm_b = (b0 + b1) / 2 - pb
    travel[0] = 1.0
    travel[1] = 1j
    travel[2] = 1j * arm_a
    travel[3] = -1.0
    travel[4] = -1j
    travel[5] = -1j * arm_b
    for p in range(6):
        twist[p] = 0j
    twist[2] = 2j * s4 / center
    twist[5] = 2j * s5 / center
    area = abs(a1 - a0) * abs(b1 - b0)
    for p in range(6):
        gradient[p, n] = area * (slope * travel[p]).real
    gradient[2, n] -= area * (1j * s2).real
    gradient[5, n] -= area * (1j * s3).real
    for p in range(6):
        for q in range(6):
            entry = bend * travel[p] * travel[q] + twist[p] * travel[q] + travel[p] * twist[q]
            hessian[p, q, n] = area * entry.real
    hessian[2, 2, n] += area * (s6 - slope * arm_a).real
    hessian[5, 5, n] += area * (s8 + slope * arm_b).real
    hessian[2, 5, n] += area * s7.real
    hessian[5, 2, n] += area * s7.real


@numba.njit(cache=True, nogil=True, fastmath=FLAGS)
def _derivatives(a0, a1, b0, b1, pa, pb, gradient, hessian):
    ap = np.empty(ORDER + 1, np.complex128)
    bp = np.empty(ORDER + 1, np.complex128)
    first = np.empty(3, np.complex128)
    second = np.empty((3, 3), np.complex128)
    pull = np.empty(6, np.complex128)
    travel = np.empty(6, np.complex128)
    twist = np.empty(6, np.complex128)
    for n in range(a0.size):
        center = (a0[n] + a1[n] - b0[n] - b1[n]) / 2
        if abs(a1[n] - a0[n]) + abs(b1[n] - b0[n]) > 2 * REACH * abs(center):
            _closed(
                a0[n],
                a1[n],
                b0[n],
                b1[n],
                center,
                pa[n],
                pb[n],
                first,
                second,
                pull,
                gradient,
                hessian,
                n,
            )
        else:
            _series(
                a0[n],
                a1[n],
                b0[n],
                b1[n],
                center,
                pa[n],
                pb[n],
                ap,
                bp,
                travel,
                twist,
                gradient,
                hessian,
                n,
            )


def _flat(shape: tuple[int, ...], *arrays: Complex) -> list[Complex]:
    return [
        np.ascontiguousarray(np.broadcast_to(array, shape), dtype=np.complex128).reshape(-1)
        for array in arrays
    ]


def pair_integrals(a0: Complex, a1: Complex, b0: Complex, b1: Complex) -> Array:
    """Double integral of ln|r - r'| over element pairs, equation (64), as in ``kernel``."""
    shape = np.broadcast_shapes(np.shape(a0), np.shape(a1), np.shape(b0), np.shape(b1))
    out = np.empty(math.prod(shape))
    _integrals(*_flat(shape, a0, a1, b0, b1), out)
    return out.reshape(shape)


def pair_derivatives(
    a0: Complex, a1: Complex, b0: Complex, b1: Complex, pivot_a: Complex, pivot_b: Complex
) -> tuple[Array, Array]:
    """Gradient and Hessian of the pair integral for rigid motion of both elements,
    equations (83) to (88), as in ``kernel``."""
    arrays = (a0, a1, b0, b1, pivot_a, pivot_b)
    shape = np.broadcast_shapes(*(np.shape(array) for array in arrays))
    count = math.prod(shape)
    gradient, hessian = np.empty((6, count)), np.empty((6, 6, count))
    _derivatives(*_flat(shape, *arrays), gradient, hessian)
    return gradient.reshape(6, *shape), hessian.reshape(6, 6, *shape)

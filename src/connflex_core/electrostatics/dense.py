"""Dense solution of the field of a mesh of boundary elements."""

from __future__ import annotations

from collections.abc import Callable
from contextlib import AbstractContextManager
from dataclasses import dataclass
from functools import partial

import numpy as np
import scipy.linalg
from numpy.typing import NDArray

from connflex_core.electrostatics.kernel import pair_derivatives, pair_integrals
from connflex_core.electrostatics.outlines import Mesh
from connflex_core.progress import Stage
from connflex_core.workers import Workers

Array = NDArray[np.float64]
Complex = NDArray[np.complex128]
Enter = Callable[..., AbstractContextManager[Stage]]

PAIRS_PER_PASS = 250_000
DERIVATIVE_PAIRS_PER_PASS = 40_000


@dataclass(frozen=True, eq=False)
class Charges:
    """Charges that unit potentials on the conductors put on the elements of a mesh.

    Lengths are in units of ``scale`` and charges per unit permittivity.

    Attributes:
        capacitance: n by n charge on the core elements of each conductor per volt on
            each conductor, equation (92). Symmetric when every element is in the core.
        charge: N by n charge on each element per volt on each conductor.
        taken: N by n charge on each element per volt on the core elements of each
            conductor alone; ``charge`` itself when every element is in the core.
        conductor: Conductor of each element.
        start: Start points of the elements as complex numbers.
        end: End points.
        length: Lengths of the elements.
        scale: Unit of length, in meters.
        origin: Point in meters that ``start`` and ``end`` are measured from.
        inverse: Solution of the matrix of equation (64) for given potentials.
        uniform: Solution of that matrix for one volt on every element.
        period: Period of the mesh in units of ``scale`` as a complex number, None
            when the mesh is not repeated.
    """

    capacitance: Array
    charge: Array
    taken: Array
    conductor: NDArray[np.intp]
    start: Complex
    end: Complex
    length: Array
    scale: float
    origin: Array
    inverse: Callable[[Array], Array]
    uniform: Array
    period: complex | None = None


def solve_charges(
    mesh: Mesh,
    conductors: int,
    enter: Enter,
    core: NDArray[np.bool_] | None = None,
    period: Array | None = None,
    workers: Workers | None = None,
) -> Charges:
    """Solve the field of a mesh for one volt on each of its conductors in turn.

    ``core`` marks the elements whose charge counts toward ``capacitance``; by default
    all do. With ``period``, a vector in meters, the mesh is one cell of a row repeated
    without end along it, equation (151). ``enter`` opens the stages assemble and solve:
    ``stage`` bound to a progress receiver, or ``quiet`` when the solution is one step of
    a larger stage. ``workers`` share the assembly by blocks of rows, which assemble
    advances once each.
    """
    count = len(mesh.conductor)
    corners = np.concatenate([mesh.start, mesh.end])
    low, high = corners.min(axis=0), corners.max(axis=0)
    origin = (low + high) / 2
    scale = 2 * float(np.hypot(*(high - low)))
    start, end = ((points - origin) / scale @ (1, 1j) for points in (mesh.start, mesh.end))
    length = np.abs(end - start)
    step = None if period is None else complex(*(np.asarray(period) / scale))

    potential = np.empty((count, count))
    rows = max(1, PAIRS_PER_PASS // count)
    blocks = [slice(first, first + rows) for first in range(0, count, rows)]

    def fill(block: slice) -> None:
        if step is None:
            integral = pair_integrals(start[block, None], end[block, None], start, end)
        else:
            from connflex_core.electrostatics.periodic import periodic_integrals

            integral = periodic_integrals(start[block, None], end[block, None], start, end, step)
        potential[block] = -integral / (2 * np.pi * length[block, None] * length)

    with enter("assemble", total=len(blocks)) as current:
        (workers or Workers(1)).run([partial(fill, block) for block in blocks], current.advance)

    with enter("solve"):
        try:
            if step is None:
                cholesky = scipy.linalg.cho_factor(potential, lower=True, overwrite_a=True)

                def inverse(right: Array) -> Array:
                    solved: Array = scipy.linalg.cho_solve(cholesky, right)
                    return solved

            else:
                pivoted = scipy.linalg.lu_factor(potential, overwrite_a=True)

                def inverse(right: Array) -> Array:
                    solved: Array = scipy.linalg.lu_solve(pivoted, right)
                    return solved

        except scipy.linalg.LinAlgError as exc:
            raise ValueError("outlines cross or elements coincide") from exc
        on = (mesh.conductor[:, None] == np.arange(conductors)).astype(float)
        uniform: Array = inverse(np.ones(count))

        def neutral(potentials: Array) -> Array:
            driven: Array = inverse(potentials)
            shift: Array = np.outer(uniform, uniform @ potentials) / float(uniform.sum())
            return driven - shift

        charge = neutral(on)
        inside = on if core is None else on * core[:, None]
        taken = charge if core is None else neutral(inside)
    return Charges(
        capacitance=inside.T @ charge,
        charge=charge,
        taken=taken,
        conductor=mesh.conductor,
        start=start,
        end=end,
        length=length,
        scale=scale,
        origin=origin,
        inverse=inverse,
        uniform=uniform,
        period=step,
    )


def charge_derivatives(
    charges: Charges, body: NDArray[np.intp], pivots: Array, enter: Enter
) -> tuple[Array, Array]:
    """Gradient and Hessian of ``charges.capacitance`` for rigid motion of bodies.

    ``body`` gives the body of each conductor and ``pivots`` the point in meters each
    body turns about. Returns n by n by m by 3 first and n by n by m by 3 by m by 3
    second derivatives per unit permittivity, in meters and radians, equations (76),
    (84), and (95). ``enter`` opens the stages derivatives and curvature.
    """
    count, conductors = charges.charge.shape
    rigid = len(pivots)
    start, end, length = charges.start, charges.end, charges.length
    whole = charges.taken is charges.charge
    right = charges.charge / length[:, None]
    left = right if whole else charges.taken / length[:, None]
    carried = [np.flatnonzero(body[charges.conductor] == k) for k in range(rigid)]
    pivot = (pivots - charges.origin) / charges.scale @ (1, 1j)
    shape = (conductors, conductors, 3 * rigid)
    gradient = np.zeros(shape)
    hessian = np.zeros((*shape, 3 * rigid))
    swing = np.zeros((count, conductors, 3 * rigid))
    counter = swing if whole else np.zeros_like(swing)
    pairs = [(k, other) for k in range(rigid) for other in range(k + 1, rigid)]
    with enter("derivatives", total=len(pairs)) as current:
        for k, other in pairs:
            mine, theirs = carried[k], carried[other]
            freedom = np.r_[3 * k : 3 * k + 3, 3 * other : 3 * other + 3]
            rows = max(1, DERIVATIVE_PAIRS_PER_PASS // max(theirs.size, 1))
            for first in range(0, mine.size if theirs.size else 0, rows):
                part = mine[first : first + rows]
                ends = (start[part, None], end[part, None], start[theirs], end[theirs])
                if charges.period is None:
                    slopes, bends = pair_derivatives(*ends, pivot[k], pivot[other])
                else:
                    from connflex_core.electrostatics.periodic import periodic_derivatives

                    slopes, bends = periodic_derivatives(
                        *ends, pivot[k], pivot[other], charges.period
                    )
                for p, slope in enumerate(slopes):
                    onto, back = slope @ right[theirs], slope.T @ right[part]
                    gradient[:, :, freedom[p]] += left[part].T @ onto
                    swing[part, :, freedom[p]] += onto / length[part, None]
                    swing[theirs, :, freedom[p]] += back / length[theirs, None]
                    if not whole:
                        gradient[:, :, freedom[p]] += left[theirs].T @ back
                        counter[part, :, freedom[p]] += slope @ left[theirs] / length[part, None]
                        counter[theirs, :, freedom[p]] += (
                            slope.T @ left[part] / length[theirs, None]
                        )
                    for q in range(6):
                        bend = bends[p, q]
                        hessian[:, :, freedom[p], freedom[q]] += left[part].T @ bend @ right[theirs]
                        if not whole:
                            hessian[:, :, freedom[p], freedom[q]] += (
                                left[theirs].T @ bend.T @ right[part]
                            )
            current.advance()
    if whole:
        gradient = gradient + gradient.transpose(1, 0, 2)
        hessian = hessian + hessian.transpose(1, 0, 2, 3)
    gradient /= 2 * np.pi
    hessian /= 2 * np.pi

    with enter("curvature"):
        flat = swing.reshape(count, -1) / (-2 * np.pi)
        settled = charges.inverse(flat)
        settled -= np.outer(charges.uniform, charges.uniform @ flat) / charges.uniform.sum()
        across = flat if whole else counter.reshape(count, -1) / (-2 * np.pi)
        through = (across.T @ settled).reshape(*shape[1:], *shape[1:]).transpose(0, 2, 1, 3)
        hessian += through + through.transpose(0, 1, 3, 2)

    unit = np.tile([1 / charges.scale, 1 / charges.scale, 1.0], rigid)
    gradient *= unit
    hessian *= unit[:, None] * unit
    return (
        gradient.reshape(conductors, conductors, rigid, 3),
        hessian.reshape(conductors, conductors, rigid, 3, rigid, 3),
    )

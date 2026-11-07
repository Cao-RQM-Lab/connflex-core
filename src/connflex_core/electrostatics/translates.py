"""Field solutions shared by translated copies of geometry: their keys and pivots."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from connflex_core.electrostatics.outlines import Mesh

Array = NDArray[np.float64]
Index = NDArray[np.intp]

KEY_GRID = 1e-11


def canonical(
    mesh: Mesh, elements: Index, body: Index, reference: Array, spots: dict[int, complex]
) -> tuple[Index, Index, Index, Index, bytes]:
    """An order of ``elements`` that translates share, their conductors and bodies by
    first appearance in it, the rank among those bodies of each conductor's body, and a
    key of that geometry relative to ``reference``, with the point ``spots`` gives a
    body, if any, as part of it."""
    spans = np.round(
        np.hstack([mesh.start[elements] - reference, mesh.end[elements] - reference]) / KEY_GRID
    ).astype(np.int64)
    order = np.lexsort(spans.T[::-1])
    owners = mesh.conductor[elements[order]]
    _, first = np.unique(owners, return_index=True)
    present = owners[np.sort(first)]
    owned = body[present]
    _, first_body = np.unique(owned, return_index=True)
    moving = owned[np.sort(first_body)]
    labels = {int(g): k for k, g in enumerate(present)}
    kinds = np.array([list(moving).index(b) for b in owned], dtype=np.intp)
    turned = np.array(
        [
            (
                np.round(
                    (np.array([spots[int(b)].real, spots[int(b)].imag]) - reference) / KEY_GRID
                ).astype(np.int64)
                if int(b) in spots
                else np.zeros(2, dtype=np.int64)
            )
            for b in moving
        ]
    )
    key = (
        spans[order].tobytes()
        + np.array([labels[int(g)] for g in owners], dtype=np.int64).tobytes()
        + kinds.astype(np.int64).tobytes()
        + turned.tobytes()
    )
    return order, present, moving, kinds, key


def taylor(value: Array, slope: Array, bend: Array, step: Array) -> tuple[Array, Array]:
    """Value and gradient a step away from a kept solution, to second order with the
    Hessian held, equation (149): ``slope`` holds the parameters in its last axis and
    ``bend`` in its last two."""
    turned = bend @ step
    return value + slope @ step + 0.5 * turned @ step, slope + turned


def repivot(slope: Array, bend: Array, arms: Array) -> tuple[Array, Array]:
    """Gradient and Hessian for rigid motion of bodies about their pivots, from those
    about points ``arms`` away from the pivots, equation (155).

    ``slope`` and ``bend`` hold travel along x and y and turn of each body in order in
    their last one and two axes; ``arms`` holds one vector per body, in meters.
    """
    rigid = len(arms)
    jacobian = np.eye(3 * rigid)
    for j, (x, y) in enumerate(arms):
        jacobian[3 * j, 3 * j + 2] = -y
        jacobian[3 * j + 1, 3 * j + 2] = x
    moved_slope: Array = slope @ jacobian
    moved_bend: Array = np.einsum("pq,...pr,rs->...qs", jacobian, bend, jacobian)
    for j, (x, y) in enumerate(arms):
        moved_bend[..., 3 * j + 2, 3 * j + 2] -= slope[..., 3 * j] * x + slope[..., 3 * j + 1] * y
    return moved_slope, moved_bend

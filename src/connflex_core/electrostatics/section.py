"""Field of a cut across a wall of the device layer: what lies above and below it."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import pairwise

import numpy as np
from numpy.typing import NDArray

from connflex_core.electrostatics.dense import charge_derivatives, solve_charges
from connflex_core.electrostatics.outlines import Mesh, edge_fractions
from connflex_core.progress import quiet

Array = NDArray[np.float64]
Index = NDArray[np.intp]

WALL_ELEMENTS = 8
HANDLE_MARGIN = 4.0


@dataclass(frozen=True, eq=False)
class Cut:
    """Conductors that a line along the outward normal of a wall crosses, in the frame
    of that wall, equation (157).

    Attributes:
        start: k starts of the intervals of silicon along the line, in meters from the
            wall, ascending; the wall's own interval ends at 0 and the gap ahead of the
            wall is positive.
        end: k ends of those intervals.
        label: Conductor of each interval: 0 for the wall's, the others numbered from 1
            by first appearance.
        body: Body of each conductor, numbered by first appearance from 0, the handle's
            last when there is one.
        handle: Whether the handle lies beneath, as the last conductor.
        ahead: Whether the cut counts what lies ahead of the wall; without, it holds the
            wall's interval and what lies behind it alone.
    """

    start: Array
    end: Array
    label: Index
    body: Index
    handle: bool
    ahead: bool

    @property
    def conductors(self) -> int:
        """Number of conductors, the handle included."""
        return int(np.max(self.label, initial=-1)) + 1 + int(self.handle)


@dataclass(frozen=True, eq=False)
class Section:
    """Charge that the field of a cut puts on the core of its wall, less what the
    layout-plane field and the plates beneath already hold, per unit length and
    permittivity, equations (158) and (159).

    Attributes:
        value: Excess charge on the core per volt on each conductor of the cut.
        gradient: Its derivatives in the travel of each body of the cut along the normal
            and upward, conductors by bodies by 2.
        hessian: Second derivatives, conductors by bodies by 2 by bodies by 2.
    """

    value: Array
    gradient: Array
    hessian: Array


def solve_section(cut: Cut, thickness: float, box: float, reach: float, grading: float) -> Section:
    """Field of the device layer along ``cut``, ``thickness`` tall, with the handle
    ``box`` beneath it when the cut has one, equations (157) to (159).

    Every interval is a block of silicon of the full thickness; the handle is a slab
    under all of them that reaches ``HANDLE_MARGIN`` times ``thickness`` past the cut
    of ``reach`` on either side. Elements are at most a ``WALL_ELEMENTS``-th of the
    thickness on walls and on faces away from features, at most the smaller of that
    and ``box`` on faces where they pass the wall or the edge of a block, and shorten
    toward corners by ``grading``. The slab is so meshed under every block edge, and
    the derivatives hold that mesh where it is.
    """
    devices = int(np.max(cut.label)) + 1
    conductors = cut.conductors
    wall = int(np.flatnonzero((cut.label == 0) & (np.abs(cut.end) <= 1e-12 * reach))[0])
    width = -float(cut.start[wall])
    upright = thickness / WALL_ELEMENTS
    flat = min(box, upright)
    rectangles = [
        (
            float(low),
            0.0,
            float(high),
            thickness,
            int(label),
            (low + high) / 2 if k == wall else None,
        )
        for k, (low, high, label) in enumerate(zip(cut.start, cut.end, cut.label, strict=True))
    ]
    if cut.handle:
        margin = reach + HANDLE_MARGIN * thickness
        rectangles.append((-margin, -2 * box, margin, -box, devices, 0.0))
    mesh = _rectangles(rectangles, flat, upright, grading)
    tolerance = 1e-9 * reach
    middle = mesh.center
    own = mesh.conductor == 0
    upright_wall = (
        own & (np.abs(mesh.start[:, 0]) <= tolerance) & (np.abs(mesh.end[:, 0]) <= tolerance)
    )
    faces = (
        own
        & (np.abs(mesh.start[:, 1] - mesh.end[:, 1]) <= tolerance)
        & (middle[:, 0] > -width / 2)
        & (middle[:, 0] < tolerance)
    )
    charges = solve_charges(mesh, conductors, quiet, upright_wall | faces)
    bodies = int(np.max(cut.body)) + 1
    # The handle is a plane: along the cut it goes with the wall, which it is centered
    # on, and only its rise is its body's own.
    moving = cut.body.copy()
    extra = int(cut.handle)
    if cut.handle:
        moving[-1] = bodies
    slope, bend = charge_derivatives(charges, moving, np.zeros((bodies + extra, 2)), quiet)
    to = np.add.outer(2 * np.arange(bodies + extra), np.arange(2))
    if cut.handle:
        to[bodies] = [2 * int(cut.body[0]), 2 * int(cut.body[-1]) + 1]
    to = to.ravel()
    raw = slope[0, :, :, :2].reshape(conductors, -1)
    raw_bend = bend[0, :, :, :2, :, :2].reshape(conductors, len(to), len(to))
    gradient: Array = np.zeros((conductors, 2 * bodies))
    hessian: Array = np.zeros((conductors, 2 * bodies, 2 * bodies))
    np.add.at(gradient, (slice(None), to), raw)
    np.add.at(hessian, (slice(None), to[:, None], to[None, :]), raw_bend)
    gradient = gradient.reshape(conductors, bodies, 2)
    hessian = hessian.reshape(conductors, bodies, 2, bodies, 2)
    value = charges.capacitance[0].copy()
    ahead = cut.start > tolerance
    if cut.ahead and ahead.any():
        first = int(np.flatnonzero(ahead)[0])
        _subtract(
            value,
            gradient,
            hessian,
            0,
            int(cut.label[first]),
            cut.body,
            0,
            thickness,
            float(cut.start[first]),
        )
    if cut.handle:
        _subtract(value, gradient, hessian, 0, devices, cut.body, 1, width / 2, box)
    return Section(value, gradient, hessian)


def _subtract(
    value: Array,
    gradient: Array,
    hessian: Array,
    near: int,
    far: int,
    body: Index,
    axis: int,
    height: float,
    gap: float,
) -> None:
    """Take ``height / gap`` between conductors ``near`` and ``far`` out of the row of
    ``near``, the gap closing as ``far`` travels toward ``near`` along ``axis``:
    equation (159): the layout plane along the normal, a plate upward."""
    share = height / gap
    value[near] -= share
    value[far] += share
    one, other = int(body[near]), int(body[far])
    if one == other:
        return
    slope, bend = height / gap**2, 2 * height / gap**3
    direction = 1.0 if axis == 0 else -1.0
    for row, sign in ((near, -1.0), (far, 1.0)):
        gradient[row, other, axis] -= sign * direction * slope
        gradient[row, one, axis] += sign * direction * slope
        hessian[row, other, axis, other, axis] += sign * bend
        hessian[row, one, axis, one, axis] += sign * bend
        hessian[row, one, axis, other, axis] -= sign * bend
        hessian[row, other, axis, one, axis] -= sign * bend


def _rectangles(
    rectangles: list[tuple[float, float, float, float, int, float | None]],
    flat: float,
    upright: float,
    grading: float,
) -> Mesh:
    """Elements on the outlines of axis-aligned rectangles, counterclockwise, at most
    ``upright`` long on vertical sides and on horizontal sides away from features, at
    most ``flat`` long on horizontal sides where they pass under or over the edge of
    another rectangle or a split point, and shortening to the smaller over ``grading``
    at corners. A rectangle with a split point gets a node there on its top and bottom.
    Blocks do not overlap, so only a slab that spans them is meshed to their edges."""
    corner = min(flat, upright) / grading
    edges = sorted({u for low_u, _, high_u, *_ in rectangles for u in (low_u, high_u)})
    starts, ends, owners = [], [], []
    for low_u, low_z, high_u, high_z, owner, split in rectangles:
        inner = {u for u in edges if low_u < u < high_u}
        if split is not None:
            inner.add(split)
        stops = [low_u, *sorted(inner), high_u]
        bottom = list(pairwise(stops))
        top = [(b, a) for a, b in reversed(bottom)]
        ends_of = (low_u, high_u)
        sides = [((a, low_z), (b, low_z), a in ends_of, b in ends_of) for a, b in bottom]
        sides.append(((high_u, low_z), (high_u, high_z), True, True))
        sides += [((a, high_z), (b, high_z), a in ends_of, b in ends_of) for a, b in top]
        sides.append(((low_u, high_z), (low_u, low_z), True, True))
        for first, last, left, right in sides:
            begin, finish = np.array(first), np.array(last)
            length = float(np.linalg.norm(finish - begin))
            fractions = edge_fractions(
                length, corner if left else flat, corner if right else flat, upright
            )
            nodes = begin + fractions[:, None] * (finish - begin)
            starts.append(nodes[:-1])
            ends.append(nodes[1:])
            owners.append(np.full(len(fractions) - 1, owner, dtype=np.intp))
    return Mesh(np.concatenate(starts), np.concatenate(ends), np.concatenate(owners))

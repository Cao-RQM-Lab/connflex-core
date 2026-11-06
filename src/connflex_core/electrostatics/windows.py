"""Field of a mesh too large for one dense solution, solved in overlapping windows."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import NamedTuple

import numpy as np
import shapely
from numpy.typing import NDArray
from shapely import STRtree

from connflex_core.electrostatics.dense import Charges, charge_derivatives, solve_charges
from connflex_core.electrostatics.outlines import Mesh
from connflex_core.electrostatics.translates import KEY_GRID, canonical, repivot
from connflex_core.progress import quiet

Array = NDArray[np.float64]
Index = NDArray[np.intp]

FIRST_OVERLAP = 4.0


def cut_cores(
    mesh: Mesh,
    most: int,
    elements: Index | None = None,
    crowd: Callable[[Index], int] | None = None,
) -> list[Index]:
    """Split elements of a mesh into compact groups of at most ``most``, equation (91).

    ``elements`` are the ones to split, all of the mesh by default. A group is halved at
    the median of its element centers along the longer side of their bounding box until
    it is small enough, and with ``crowd``, which counts the elements of the first window
    around a group, until that window holds at most twice ``most``. Centers are ordered
    as rounded for keys of translates, so that groups that are translates are cut alike.
    """
    center = (mesh.start + mesh.end) / 2
    every = np.arange(len(mesh.conductor)) if elements is None else np.asarray(elements)
    cores, pending = [], [every]
    while pending:
        group = pending.pop()
        crowded = crowd is not None and len(group) > 1 and crowd(group) > 2 * most
        if len(group) <= most and not crowded:
            cores.append(group)
            continue
        points = center[group]
        grid = np.round((points - points.min(axis=0)) / KEY_GRID).astype(np.int64)
        axis = int(np.argmax(np.ptp(grid, axis=0)))
        order = group[np.lexsort((grid[:, 1 - axis], grid[:, axis]))]
        pending += [order[len(group) // 2 :], order[: len(group) // 2]]
    return cores


def smallest_window(mesh: Mesh, core: Index, size: float, reach: float | None) -> int:
    """Number of elements in the first window that ``solve_windows`` tries around ``core``."""
    center = (mesh.start + mesh.end) / 2
    low, high = center[core].min(axis=0), center[core].max(axis=0)
    overlap = FIRST_OVERLAP * size if reach is None else min(FIRST_OVERLAP * size, reach)
    inside = np.all((center >= low - overlap) & (center <= high + overlap), axis=1)
    return int(inside.sum())


class _Window(NamedTuple):
    """Elements of a window in an order that translates share, its conductors and bodies
    by first appearance in that order, and which elements are its core."""

    elements: Index
    present: Index
    moving: Index
    mine: NDArray[np.bool_]


@dataclass
class _Solved:
    """Charge on the core of a window per volt on each of its conductors, in an order
    that translates share, and its derivatives about the corner of the core."""

    capacitance: Array
    charge: Array
    derivatives: tuple[Array, Array] | None = None


def solve_windows(
    mesh: Mesh,
    body: Index,
    pivots: Array,
    cores: list[Index],
    size: float,
    reach: float | None,
    tolerance: float,
    limit: int,
    advance: Callable[[], None],
    windows: list[Index] | None = None,
    pairs: Index | None = None,
) -> tuple[Array, Array, Array, Array, list[Index]]:
    """Sum over windows of the charge that each solution puts on its core.

    A window holds a core and every element within an overlap of its bounding box,
    equations (91) to (96). The overlap starts at ``FIRST_OVERLAP`` times ``size`` and
    doubles until the charge on the core changes by less than ``tolerance``, the overlap
    reaches ``reach``, or the next window would hold more than ``limit`` elements.
    With ``windows``, the elements of each window as an earlier call returned them, no
    overlap is searched: a mesh whose bodies have moved keeps the windows it had. With
    ``pairs``, pairs of conductors as rows, a window is solved only when it holds both
    conductors of a pair and its core lies on one of them; otherwise when it holds two
    conductors. A window that is a translate of one solved before, with its core,
    conductors, and bodies alike, takes that solution, its derivatives moved to its own
    pivots by equation (155).
    Returns capacitance, gradient, Hessian, and element charges per unit permittivity,
    in the layout of ``SidewallCapacitance`` and before capacitance is made symmetric,
    and the elements of each window, none for a window that was not solved.
    ``advance`` is called once per core.
    """
    count, conductors, rigid = len(mesh.conductor), len(body), len(pivots)
    center = (mesh.start + mesh.end) / 2
    tree = STRtree(np.asarray(shapely.points(center)))
    capacitance = np.zeros((conductors, conductors))
    gradient = np.zeros((conductors, conductors, rigid, 3))
    hessian = np.zeros((conductors, conductors, rigid, 3, rigid, 3))
    charge = np.zeros((count, conductors))
    three = np.arange(3)
    chosen: list[Index] = []
    known: dict[bytes, _Solved] = {}
    rank = np.zeros(conductors, dtype=np.intp)

    def useful(present: Index, core: Index) -> bool:
        if pairs is None:
            return len(present) > 1
        on = np.isin(pairs, np.unique(mesh.conductor[core])).any(axis=1)
        return bool(np.any(np.isin(pairs, present).all(axis=1) & on))

    def charges_of(elements: Index, present: Index, mine: NDArray[np.bool_]) -> Charges:
        rank[present] = np.arange(len(present))
        part = Mesh(mesh.start[elements], mesh.end[elements], rank[mesh.conductor[elements]])
        return solve_charges(part, len(present), quiet, mine)

    def solve(inside: Index, core: Index) -> tuple[_Window, _Solved, Charges | None]:
        corner = center[core].min(axis=0)
        order, present, moving, key = canonical(mesh, inside, body, corner, {})
        elements = inside[order]
        mine = np.isin(elements, core)
        key += mine.tobytes()
        charges = None
        if key not in known:
            charges = charges_of(elements, present, mine)
            known[key] = _Solved(charges.capacitance, charges.charge[mine])
        return _Window(elements, present, moving, mine), known[key], charges

    for number, core in enumerate(cores):
        low, high = center[core].min(axis=0), center[core].max(axis=0)
        overlap = FIRST_OVERLAP * size if reach is None else min(FIRST_OVERLAP * size, reach)
        held: tuple[_Window, _Solved, Charges | None] | None = None
        share = np.zeros((conductors, conductors))
        smallest = True
        if windows is not None and len(windows[number]):
            held = solve(windows[number], core)
        while windows is None:
            inside = np.sort(tree.query(shapely.box(*(low - overlap), *(high + overlap))))
            if len(inside) > limit:
                if smallest:
                    raise ValueError(
                        f"a window needs {len(inside)} boundary elements, more than "
                        f"max_elements = {limit}: raise max_elements or element_size, or "
                        f"lower window_elements"
                    )
                break
            smallest = False
            present = np.unique(mesh.conductor[inside])
            before, share = share, np.zeros((conductors, conductors))
            settled = False
            if useful(present, core):
                found = solve(inside, core)
                ordered = found[0].present
                share[np.ix_(ordered, ordered)] = found[1].capacitance
                settled = held is not None and bool(
                    np.max(np.abs(share - before)) < tolerance * np.max(np.abs(share))
                )
                held = found
            if settled or len(inside) == count or (reach is not None and overlap >= reach):
                break
            overlap = 2 * overlap if reach is None else min(2 * overlap, reach)
        if held is not None:
            (elements, present, moving, mine), solved, charges = held
            corner, here, turning = center[core].min(axis=0), len(present), 3 * len(moving)
            if solved.derivatives is None:
                charges = charges or charges_of(elements, present, mine)
                place = {int(b): j for j, b in enumerate(moving)}
                kinds = np.array([place[int(b)] for b in body[present]], dtype=np.intp)
                slope, bend = charge_derivatives(
                    charges, kinds, np.tile(corner, (len(moving), 1)), quiet
                )
                solved.derivatives = (
                    slope.reshape(here, here, turning),
                    bend.reshape(here, here, turning, turning),
                )
            slope, bend = repivot(*solved.derivatives, corner - pivots[moving])
            capacitance[np.ix_(present, present)] += solved.capacitance
            gradient[np.ix_(present, present, moving)] += slope.reshape(here, here, -1, 3)
            hessian[np.ix_(present, present, moving, three, moving, three)] += bend.reshape(
                here, here, len(moving), 3, len(moving), 3
            )
            charge[np.ix_(elements[mine], present)] = solved.charge
        chosen.append(np.zeros(0, dtype=np.intp) if held is None else held[0].elements)
        advance()
    return capacitance, gradient, hessian, charge, chosen

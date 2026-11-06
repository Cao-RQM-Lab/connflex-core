"""Field of a mesh too large for one dense solution, solved in overlapping windows."""

from __future__ import annotations

import threading
from collections.abc import Callable, Hashable, Sequence
from dataclasses import dataclass, field
from functools import partial
from typing import NamedTuple

import numpy as np
import shapely
from numpy.typing import NDArray
from shapely import STRtree

from connflex_core.electrostatics.dense import Charges, charge_derivatives, solve_charges
from connflex_core.electrostatics.outlines import Mesh
from connflex_core.electrostatics.translates import KEY_GRID, canonical, repivot
from connflex_core.progress import quiet
from connflex_core.workers import Memo, Workers

Array = NDArray[np.float64]
Index = NDArray[np.intp]

FIRST_OVERLAP = 4.0
TRAVEL_MARGIN = 0.25


@dataclass(eq=False)
class WindowSets:
    """Elements of each window wherever it was chosen, equation (164), by key of its
    core.

    Attributes:
        overlap: Overlap of each window, as found at rest or as first given.
        sets: Sets of elements chosen for each window so far, the one at rest first,
            each with the overlap it was chosen with; empty for a window not solved.
        last: Position in ``sets`` of the set each window took last.
    """

    overlap: dict[Hashable, float] = field(default_factory=dict)
    sets: dict[Hashable, list[tuple[Index, float]]] = field(default_factory=dict)
    last: dict[Hashable, int] = field(default_factory=dict)


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


def center_tree(mesh: Mesh) -> STRtree:
    """Index of element centers, for the elements within a box."""
    return STRtree(np.asarray(shapely.points((mesh.start + mesh.end) / 2)))


def smallest_window(
    mesh: Mesh, tree: STRtree, core: Index, size: float, reach: float | None
) -> int:
    """Number of elements in the first window that ``solve_windows`` tries around
    ``core``; ``tree`` indexes the centers of ``mesh`` as ``center_tree`` builds it."""
    center = (mesh.start[core] + mesh.end[core]) / 2
    low, high = center.min(axis=0), center.max(axis=0)
    overlap = FIRST_OVERLAP * size if reach is None else min(FIRST_OVERLAP * size, reach)
    return len(tree.query(shapely.box(*(low - overlap), *(high + overlap))))


class _Window(NamedTuple):
    """Elements of a window in an order that translates share, its conductors and bodies
    by first appearance in that order, and which elements are its core."""

    elements: Index
    present: Index
    moving: Index
    mine: NDArray[np.bool_]


@dataclass(eq=False)
class _Solved:
    """Charge on the core of a window per volt on each of its conductors, in an order
    that translates share, and its derivatives about the corner of the core, found once
    under ``lock``; ``charges`` keeps the solution for them, in the coordinates of the
    window that solved it, whose core corner is ``corner``, until they are found or that
    window moves on to a wider one."""

    capacitance: Array
    charge: Array
    charges: Charges | None
    corner: Array
    derivatives: tuple[Array, Array] | None = None
    lock: threading.Lock = field(default_factory=threading.Lock)


class _Found(NamedTuple):
    """What the window of one core gave: its overlap, its elements, and when it was
    solved its conductors, its bodies, its derivatives at the pivots, its capacitance,
    its core elements, and their charges."""

    overlap: float
    elements: Index
    solved: tuple[Index, Index, Array, Array, Array, Index, Array] | None


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
    windows: WindowSets | None = None,
    pairs: Index | None = None,
    keys: Sequence[Hashable] | None = None,
    workers: Workers | None = None,
) -> tuple[Array, Array, Array, Array, WindowSets]:
    """Sum over windows of the charge that each solution puts on its core.

    A window holds a core and every element within an overlap of its bounding box,
    equations (91) to (96). The overlap starts at ``FIRST_OVERLAP`` times ``size`` and
    doubles until the charge on the core changes by less than ``tolerance``, the overlap
    reaches ``reach``, or the next window would hold more than ``limit`` elements.
    With ``windows``, as an earlier call returned them, no overlap is searched: a mesh
    whose bodies have moved keeps a set of elements for each window that holds what lies
    within the overlap less ``TRAVEL_MARGIN`` of it, and otherwise takes the elements
    within the overlap where they are now, equation (164). ``keys`` name the cores, by
    default their positions; a core new to ``windows`` takes ``reach`` as its overlap.
    With ``pairs``, pairs of conductors as rows, a window is solved only when it holds
    both conductors of a pair and its core lies on one of them; otherwise when it holds
    two conductors. A window that is a translate of one solved before, with its core,
    conductors, and bodies alike, takes that solution, its derivatives moved to its own
    pivots by equation (155). ``workers`` run the windows, one piece each; sums are
    taken in the order of ``cores``.
    Returns capacitance, gradient, Hessian, and element charges per unit permittivity,
    in the layout of ``SidewallCapacitance`` and before capacitance is made symmetric,
    and the sets of elements of the windows. ``advance`` is called once per core.
    """
    count, conductors, rigid = len(mesh.conductor), len(body), len(pivots)
    center = (mesh.start + mesh.end) / 2
    tree = center_tree(mesh)
    capacitance = np.zeros((conductors, conductors))
    gradient = np.zeros((conductors, conductors, rigid, 3))
    hessian = np.zeros((conductors, conductors, rigid, 3, rigid, 3))
    charge = np.zeros((count, conductors))
    three = np.arange(3)
    known: Memo[bytes, _Solved] = Memo()

    def useful(present: Index, core: Index) -> bool:
        if pairs is None:
            return len(present) > 1
        on = np.isin(pairs, np.unique(mesh.conductor[core])).any(axis=1)
        return bool(np.any(np.isin(pairs, present).all(axis=1) & on))

    def around(core: Index, overlap: float) -> Index:
        low, high = center[core].min(axis=0), center[core].max(axis=0)
        return np.sort(tree.query(shapely.box(*(low - overlap), *(high + overlap))))

    def follow(key: Hashable, core: Index, sets: WindowSets) -> Index:
        """Set of elements that window ``key`` takes where the bodies are now."""
        options = sets.sets.setdefault(key, [])

        def holds(elements: Index, overlap: float) -> bool:
            needed = around(core, (1 - TRAVEL_MARGIN) * overlap)
            if len(elements):
                return bool(np.isin(needed, elements).all())
            return not useful(np.unique(mesh.conductor[needed]), core)

        last = sets.last.get(key, 0)
        for k in sorted(range(len(options)), key=lambda k: (k != last, -k)):
            if holds(*options[k]):
                sets.last[key] = k
                return options[k][0]
        overlap = sets.overlap.setdefault(key, FIRST_OVERLAP * size if reach is None else reach)
        inside = around(core, overlap)
        while len(inside) > limit and overlap > FIRST_OVERLAP * size:
            overlap /= 2
            inside = around(core, overlap)
        if len(inside) > limit:
            raise ValueError(
                f"a window needs {len(inside)} boundary elements where the bodies are, more "
                f"than max_elements = {limit}: raise max_elements or element_size, or lower "
                f"window_elements"
            )
        if not useful(np.unique(mesh.conductor[inside]), core):
            inside = inside[:0]
        options.append((inside, overlap))
        sets.last[key] = len(options) - 1
        return inside

    def charges_of(elements: Index, present: Index, mine: NDArray[np.bool_]) -> Charges:
        rank = np.zeros(conductors, dtype=np.intp)
        rank[present] = np.arange(len(present))
        part = Mesh(mesh.start[elements], mesh.end[elements], rank[mesh.conductor[elements]])
        return solve_charges(part, len(present), quiet, mine)

    def solve(inside: Index, core: Index) -> tuple[_Window, _Solved]:
        corner = center[core].min(axis=0)
        order, present, moving, key = canonical(mesh, inside, body, corner, {})
        elements = inside[order]
        mine = np.isin(elements, core)
        key += mine.tobytes()

        def compute() -> _Solved:
            charges = charges_of(elements, present, mine)
            return _Solved(charges.capacitance, charges.charge[mine], charges, corner)

        return _Window(elements, present, moving, mine), known.get(key, compute)

    def window_of(number: int, core: Index) -> _Found:
        overlap = FIRST_OVERLAP * size if reach is None else min(FIRST_OVERLAP * size, reach)
        held: tuple[_Window, _Solved] | None = None
        share = np.zeros((conductors, conductors))
        smallest, taken = True, 0.0
        if windows is not None:
            elements = follow(number if keys is None else keys[number], core, windows)
            if len(elements):
                held = solve(elements, core)
        while windows is None:
            inside = around(core, overlap)
            if len(inside) > limit:
                if smallest:
                    raise ValueError(
                        f"a window needs {len(inside)} boundary elements, more than "
                        f"max_elements = {limit}: raise max_elements or element_size, or "
                        f"lower window_elements"
                    )
                break
            smallest = False
            if held is None:
                taken = overlap
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
                if held is not None and held[1] is not found[1]:
                    held[1].charges = None
                held = found
                taken = overlap
            if settled or len(inside) == count or (reach is not None and overlap >= reach):
                break
            overlap = 2 * overlap if reach is None else min(2 * overlap, reach)
        if held is None:
            return _Found(taken, np.zeros(0, dtype=np.intp), None)
        (elements, present, moving, mine), solved = held
        corner, here, turning = center[core].min(axis=0), len(present), 3 * len(moving)
        with solved.lock:
            if solved.derivatives is None:
                charges, about = solved.charges, solved.corner
                if charges is None:
                    charges, about = charges_of(elements, present, mine), corner
                place = {int(b): j for j, b in enumerate(moving)}
                kinds = np.array([place[int(b)] for b in body[present]], dtype=np.intp)
                slope, bend = charge_derivatives(
                    charges, kinds, np.tile(about, (len(moving), 1)), quiet
                )
                solved.derivatives = (
                    slope.reshape(here, here, turning),
                    bend.reshape(here, here, turning, turning),
                )
                solved.charges = None
        slope, bend = repivot(*solved.derivatives, corner - pivots[moving])
        return _Found(
            taken,
            elements,
            (present, moving, slope, bend, solved.capacitance, elements[mine], solved.charge),
        )

    found = (workers or Workers(1)).run(
        [partial(window_of, number, core) for number, core in enumerate(cores)], advance
    )
    for piece in found:
        if piece.solved is None:
            continue
        present, moving, slope, bend, value, core, held = piece.solved
        here = len(present)
        capacitance[np.ix_(present, present)] += value
        gradient[np.ix_(present, present, moving)] += slope.reshape(here, here, -1, 3)
        hessian[np.ix_(present, present, moving, three, moving, three)] += bend.reshape(
            here, here, len(moving), 3, len(moving), 3
        )
        charge[np.ix_(core, present)] = held
    if windows is None:
        names = list(range(len(cores))) if keys is None else list(keys)
        windows = WindowSets(
            {key: piece.overlap for key, piece in zip(names, found, strict=True)},
            {
                key: [(piece.elements, piece.overlap)]
                for key, piece in zip(names, found, strict=True)
            },
        )
    return capacitance, gradient, hessian, charge, windows

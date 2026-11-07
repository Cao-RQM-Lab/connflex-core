"""Rows of identical cells: the periodic kernel, arrays in a mesh, and their capacitance."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field, replace

import networkx as nx
import numpy as np
import shapely
from numpy.typing import NDArray
from shapely import MultiPolygon, STRtree

from connflex_core.electrostatics.outlines import Mesh
from connflex_core.electrostatics.pairs import pair_derivatives, pair_integrals
from connflex_core.electrostatics.translates import KEY_GRID, canonical, repivot

Array = NDArray[np.float64]
Complex = NDArray[np.complex128]
Index = NDArray[np.intp]

LONGEST_MOTIF = 16
GRID = 1e-10

IMAGES = 2
_NODES, _WEIGHTS = np.polynomial.legendre.leggauss(4)
_SPOTS = (_NODES + 1) / 2
_SHARES = _WEIGHTS / 2
_SMALL = 0.25
_FAR = 20.0


def periodic_integrals(
    a0: Complex, a1: Complex, b0: Complex, b1: Complex, period: complex
) -> Array:
    """Double integral over element pairs of ln|2 sin(pi w / p)| less a constant,
    equation (151), for a row repeated without end along ``period``.

    The second element of each pair is first taken to its copy nearest the first. The
    nearest ``IMAGES`` copies on either side of that enter through the exact pair
    integral, the rest through the remainder of (152) by Gauss-Legendre quadrature on
    both elements. A constant per pair is dropped, which a neutral row does not feel.
    """
    a0, a1, b0, b1 = np.broadcast_arrays(a0, a1, b0, b1)
    shift = _nearest(a0, a1, b0, b1, period)
    b0, b1 = b0 + shift, b1 + shift
    total = np.zeros(a0.shape)
    for k in range(-IMAGES, IMAGES + 1):
        total += pair_integrals(a0, a1, b0 + k * period, b1 + k * period)
    along, reach = np.conj(period) / abs(period), abs(period)
    left, right = _spots(a0, a1), _spots(b0, b1)
    value = _remainder(along * (left[..., :, None] - right[..., None, :]), reach)[0]
    weight = np.abs(a1 - a0) * np.abs(b1 - b0)
    total += weight * np.einsum("...ab,a,b->...", value, _SHARES, _SHARES)
    return total


def periodic_derivatives(
    a0: Complex,
    a1: Complex,
    b0: Complex,
    b1: Complex,
    pivot_a: complex,
    pivot_b: complex,
    period: complex,
) -> tuple[Array, Array]:
    """Gradient and Hessian of ``periodic_integrals`` for rigid motion of both elements.

    Parameters are as for ``pair_derivatives``. Every copy of an element moves with it
    and turns about the copy of its pivot, so the row stays periodic.
    """
    a0, a1, b0, b1 = np.broadcast_arrays(a0, a1, b0, b1)
    shift = _nearest(a0, a1, b0, b1, period)
    b0, b1 = b0 + shift, b1 + shift
    copied: Complex = pivot_b + shift
    gradient = np.zeros((6, *a0.shape))
    hessian = np.zeros((6, 6, *a0.shape))
    for k in range(-IMAGES, IMAGES + 1):
        slope, bend = pair_derivatives(
            a0,
            a1,
            b0 + k * period,
            b1 + k * period,
            np.asarray(pivot_a, dtype=np.complex128),
            copied + k * period,
        )
        gradient += slope
        hessian += bend
    along, reach = np.conj(period) / abs(period), abs(period)
    first = _spots(a0, a1)[..., :, None]
    second = _spots(b0, b1)[..., None, :]
    _, rise, bow = _remainder(along * (first - second), reach)
    one = np.ones(first.shape, dtype=np.complex128) * np.ones(second.shape)
    moves = [
        along * one,
        along * 1j * one,
        along * 1j * (first - pivot_a) * np.ones(second.shape),
        -along * one,
        -along * 1j * one,
        -along * 1j * (second - copied[..., None, None]) * np.ones(first.shape),
    ]
    turns = {
        (2, 2): -along * (first - pivot_a) * np.ones(second.shape),
        (5, 5): along * (second - copied[..., None, None]) * np.ones(first.shape),
    }
    weight = np.abs(a1 - a0) * np.abs(b1 - b0)
    for p in range(6):
        gradient[p] += weight * _integrate((rise * moves[p]).real)
        for q in range(p, 6):
            term = bow * moves[p] * moves[q]
            if (p, q) in turns:
                term = term + rise * turns[p, q]
            hessian[p, q] += weight * _integrate(term.real)
            if q != p:
                hessian[q, p] = hessian[p, q]
    return gradient, hessian


def _nearest(a0: Complex, a1: Complex, b0: Complex, b1: Complex, period: complex) -> Complex:
    """Whole periods that take the second element of each pair to its copy nearest the
    first: the remainder of (152) is smooth only within a few periods."""
    along = np.conj(period) / abs(period)
    apart = (along * ((a0 + a1) - (b0 + b1)) / 2).real
    shift: Complex = np.round(apart / abs(period)) * period
    return shift


def _spots(start: Complex, end: Complex) -> Complex:
    spots: Complex = start[..., None] + (end - start)[..., None] * _SPOTS
    return spots


def _integrate(values: Array) -> Array:
    integral: Array = np.einsum("...ab,a,b->...", values, _SHARES, _SHARES)
    return integral


def _remainder(w: Complex, reach: float) -> tuple[Array, Complex, Complex]:
    """Value, slope, and curvature of the remainder of equation (152) at ``w``, the
    separation of two points turned so that the period lies along the real axis.

    Below ``_SMALL`` the terms take their series in x = pi w / p; beyond ``_FAR`` across
    the row, where the sine would overflow, they take their asymptotic form, which is
    exact to rounding there: ln|sin x| is |Im x| - ln 2, cot x is -i sign(Im x), and
    1 / sin^2 x vanishes.
    """
    x = np.pi * w / reach
    small = np.abs(x) < _SMALL
    far = np.abs(x.imag) > _FAR
    safe = np.where(small | far, 1.0, x)
    x2 = np.where(small, x * x, 0.0)
    across = np.where(far, x, 1.0)
    log_sinc = np.where(
        small,
        -x2 / 6 - x2**2 / 180 - x2**3 / 2835 - x2**4 / 37800,
        np.where(
            far,
            np.abs(across.imag) - np.log(2) - np.log(np.abs(across)),
            np.log(np.sin(safe) / safe),
        ),
    )
    cot_less = np.where(
        small,
        -x / 3 - x * x2 / 45 - 2 * x * x2**2 / 945 - x * x2**3 / 4725,
        np.where(
            far, -1j * np.sign(across.imag) - 1 / across, np.cos(safe) / np.sin(safe) - 1 / safe
        ),
    )
    csc_less = np.where(
        small,
        -1 / 3 - x2 / 15 - 2 * x2**2 / 189 - x2**3 / 675,
        np.where(far, 1 / across**2, 1 / safe**2 - 1 / np.sin(safe) ** 2),
    )
    value: Array = np.asarray(log_sinc.real, dtype=np.float64)
    slope: Complex = np.asarray(np.pi / reach * cot_less, dtype=np.complex128)
    curve: Complex = np.asarray((np.pi / reach) ** 2 * csc_less, dtype=np.complex128)
    for k in range(1, IMAGES + 1):
        square = (k * reach) ** 2
        value = value - np.log(np.abs(1 - w * w / square))
        slope = slope + 2 * w / (square - w * w)
        curve = curve + 2 * (square + w * w) / (square - w * w) ** 2
    return value, slope, curve


@dataclass(frozen=True, eq=False)
class Row:
    """Interior cells of a row of identical cells in a mesh, equation (153).

    Attributes:
        period: Vector from one cell to the next, in meters.
        cells: Elements of each interior cell, in order along ``period``.
        centers: Center of each interior cell at rest, in meters.
        middle: Position in ``cells`` of the cell solved for all of them.
        owners: Conductors whose runs make the row.
        ends: Points at rest where the runs of each owner begin and end along the row,
            owners by 2 by 2, in meters.
    """

    period: Array
    cells: tuple[Index, ...]
    centers: Array
    middle: int
    owners: Index = field(default_factory=lambda: np.zeros(0, dtype=np.intp))
    ends: Array = field(default_factory=lambda: np.zeros((0, 2, 2)))

    @property
    def elements(self) -> Index:
        """Elements of every interior cell."""
        joined: Index = np.concatenate(self.cells)
        return joined


def find_rows(regions: Sequence[MultiPolygon], mesh: Mesh, reach: float, fewest: int) -> list[Row]:
    """Rows of at least ``fewest`` identical cells whose elements are translates of each
    other, with ``reach`` of geometry beyond either end left to the windows, equation
    (153).

    A row is found where an outline repeats its edges, translated, along a run of
    corners, and its cells hold every element in the band across the row that such runs
    of two or more conductors span.
    """
    runs: dict[tuple[int, int], list[tuple[int, Array]]] = {}
    for conductor, region in enumerate(regions):
        for polygon in region.geoms:
            for ring in [polygon.exterior, *polygon.interiors]:
                for shift, points in _runs(np.asarray(ring.coords)[:-1], fewest):
                    if shift[0] < 0 or (shift[0] == 0 and shift[1] < 0):
                        shift = -shift
                    key = (round(shift[0] / GRID), round(shift[1] / GRID))
                    runs.setdefault(key, []).append((conductor, points))
    center = (mesh.start + mesh.end) / 2
    tree = STRtree(np.asarray(shapely.points(center)))
    arrays: list[Row] = []
    for key, found in runs.items():
        period = np.array(key, dtype=float) * GRID
        along = period / np.linalg.norm(period)
        across = np.array([-along[1], along[0]])
        spans = [(conductor, points @ along, points @ across) for conductor, points in found]
        linked: nx.Graph[int] = nx.Graph()
        linked.add_nodes_from(range(len(spans)))
        for i, (_, s_i, u_i) in enumerate(spans):
            for j in range(i + 1, len(spans)):
                _, s_j, u_j = spans[j]
                if (
                    np.min(s_i) < np.max(s_j)
                    and np.min(s_j) < np.max(s_i)
                    and np.min(u_i) <= np.max(u_j) + reach
                    and np.min(u_j) <= np.max(u_i) + reach
                ):
                    linked.add_edge(i, j)
        for group in nx.connected_components(linked):
            members = [spans[k] for k in group]
            owners = {conductor for conductor, _, _ in members}
            if len(owners) < 2:
                continue
            low = max(min(float(np.min(s)) for c, s, _ in members if c == o) for o in owners)
            high = min(max(float(np.max(s)) for c, s, _ in members if c == o) for o in owners)
            band = (
                min(float(np.min(u)) for *_, u in members),
                max(float(np.max(u)) for *_, u in members),
            )
            found_array = _cells(
                center,
                tree,
                along,
                across,
                period,
                members[0][1][0],
                low + reach,
                high - reach,
                band,
                mesh,
                fewest,
            )
            if found_array is not None:
                owned = sorted(owners)
                ends = []
                for owner in owned:
                    points = np.vstack([found[k][1] for k in group if found[k][0] == owner])
                    s = points @ along
                    ends.append([points[np.argmin(s)], points[np.argmax(s)]])
                arrays.append(
                    replace(found_array, owners=np.array(owned, dtype=np.intp), ends=np.array(ends))
                )
    return arrays


def _runs(points: Array, fewest: int) -> list[tuple[Array, Array]]:
    """Translations under which a ring repeats its edges over at least ``fewest``
    motifs, each with the corners of the run."""
    edges = np.round((np.roll(points, -1, axis=0) - points) / GRID).astype(np.int64)
    count = len(edges)
    found = []
    covered = np.zeros(count, dtype=bool)
    for motif in range(1, min(LONGEST_MOTIF, count // 2) + 1):
        same = np.all(edges == np.roll(edges, -motif, axis=0), axis=1)
        if same.all() or not same.any():
            continue
        start = int(np.flatnonzero(~same)[0]) + 1
        order = (np.arange(count) + start) % count
        flags = same[order]
        edges_at = np.flatnonzero(np.diff(np.concatenate([[0], flags.astype(int), [0]])))
        for first, last in zip(edges_at[::2], edges_at[1::2], strict=True):
            span = order[np.arange(first, last) % count]
            if last - first >= fewest * motif and covered[span].mean() < 0.5:
                corners = order[np.arange(first, last + motif + 1) % count]
                found.append((points[corners[motif]] - points[corners[0]], points[corners]))
                covered[span] = True
    return found


def _cells(
    center: Array,
    tree: STRtree,
    along: Array,
    across: Array,
    period: Array,
    phase: float,
    low: float,
    high: float,
    band: tuple[float, float],
    mesh: Mesh,
    fewest: int,
) -> Row | None:
    """Interior cells between ``low`` and ``high`` along the row, or None when fewer
    than ``fewest`` of them are translates of each other."""
    pitch = float(np.linalg.norm(period))
    tolerance = 1e-3 * float(np.min(np.linalg.norm(mesh.end - mesh.start, axis=1)))
    strip = tree.query(
        shapely.Polygon(
            [
                low * along + (band[0] - tolerance) * across,
                high * along + (band[0] - tolerance) * across,
                high * along + (band[1] + tolerance) * across,
                low * along + (band[1] + tolerance) * across,
            ]
        )
    )
    if not len(strip):
        return None
    turns = np.sort(((center[strip] @ along - phase) / pitch) % 1.0)
    gaps = np.diff(np.concatenate([turns, [turns[0] + 1.0]]))
    widest = int(np.argmax(gaps))
    phase += (turns[widest] + gaps[widest] / 2) % 1.0 * pitch
    first = int(np.ceil((low - phase) / pitch))
    last = int(np.floor((high - phase) / pitch))
    if last - first < fewest:
        return None
    corners = [
        low_s * along + band_u * across
        for low_s in (phase + first * pitch, phase + last * pitch)
        for band_u in (band[0] - tolerance, band[1] + tolerance)
    ]
    near = tree.query(shapely.Polygon(np.array(corners)[[0, 1, 3, 2]]).buffer(tolerance))
    s = center[near] @ along
    u = center[near] @ across
    inside = (u >= band[0] - tolerance) & (u <= band[1] + tolerance)
    position = np.floor((s - phase) / pitch).astype(int)
    cells = []
    for k in range(first, last):
        members = near[inside & (position == k)]
        cells.append(
            members[
                np.lexsort(  # positions rounded to the tolerance, so that a tie between
                    (  # mirror-image elements cannot flip from one cell to the next
                        np.round((center[members] - k * period) @ across / tolerance),
                        np.round((center[members] @ along - k * pitch) / tolerance),
                        mesh.conductor[members],
                    )
                )
            ]
        )
    reference = cells[0]
    good = []
    for offset, members in enumerate(cells):
        same = len(members) == len(reference) and np.array_equal(
            mesh.conductor[members], mesh.conductor[reference]
        )
        if same:
            moved = offset * period
            same = bool(
                np.max(np.abs(mesh.start[members] - moved - mesh.start[reference]), initial=0)
                <= tolerance
                and np.max(np.abs(mesh.end[members] - moved - mesh.end[reference]), initial=0)
                <= tolerance
            )
        good.append(same)
    runs = np.flatnonzero(np.diff(np.concatenate([[0], np.array(good, int), [0]])))
    if not len(runs):
        return None
    lengths = runs[1::2] - runs[::2]
    best = int(np.argmax(lengths))
    if lengths[best] < fewest or len(reference) == 0:
        return None
    chosen = cells[runs[2 * best] : runs[2 * best + 1]]
    centers = np.array(
        [
            phase * along
            + (first + runs[2 * best] + k + 0.5) * pitch * along
            + 0.5 * (band[0] + band[1]) * across
            for k in range(len(chosen))
        ]
    )
    return Row(period=period, cells=tuple(chosen), centers=centers, middle=len(chosen) // 2)


def near_rows(rows: Sequence[Row], mesh: Mesh, elements: Index, reach: float) -> list[Index]:
    """``elements`` in groups that translated rows with translated surroundings cut
    alike, equation (156).

    An element within twice ``reach`` and one period of the cells of a row joins the
    nearest such row. Those before its first cell form one group, those beyond its last
    another, and those beside its cells groups as long as the fewest whole cells that
    span ``reach``, counted from the first. Elements near no row form the last group.
    """
    center = (mesh.start + mesh.end) / 2
    points = center[elements]
    tree = STRtree(np.asarray(shapely.points(points)))
    nearest = np.full(len(elements), np.inf)
    owner = np.full(len(elements), -1)
    frames = []
    for k, row in enumerate(rows):
        pitch = float(np.linalg.norm(row.period))
        along = row.period / pitch
        across = np.array([-along[1], along[0]])
        low = float(row.centers[0] @ along) - pitch / 2
        high = float(row.centers[-1] @ along) + pitch / 2
        side = center[row.elements] @ across
        bottom, top = float(side.min()), float(side.max())
        frames.append((along, low, high, np.ceil(reach / pitch) * pitch))
        outline = shapely.Polygon(
            [
                s * along + u * across
                for s, u in ((low, bottom), (high, bottom), (high, top), (low, top))
            ]
        )
        near = tree.query(outline, predicate="dwithin", distance=2 * reach + pitch)
        s, u = points[near] @ along, points[near] @ across
        apart = np.hypot(
            np.maximum(np.maximum(low - s, s - high), 0.0),
            np.maximum(np.maximum(bottom - u, u - top), 0.0),
        )
        closer = apart < nearest[near]
        nearest[near[closer]] = apart[closer]
        owner[near[closer]] = k
    groups = []
    for k, (along, low, high, block) in enumerate(frames):
        mine = elements[owner == k]
        s = center[mine] @ along
        zone = np.where(s < low, -1, np.where(s > high, -2, np.floor((s - low) / block)))
        groups += [mine[zone == z] for z in np.unique(zone)]
    rest = elements[owner < 0]
    return [*groups, rest] if len(rest) else groups


def interior_cells(
    rows: Sequence[Row], mesh: Mesh, rest: Mesh, body: Index, reach: float
) -> list[NDArray[np.bool_]]:
    """Which cells of each row lie, where the bodies are now, within the overlap of its
    runs less ``reach`` at both ends, equation (165); ``mesh``, ``rest``, and ``body`` as
    for ``solve_rows``."""
    turn = _turns_of(mesh, rest, body)
    center = (mesh.start + mesh.end) / 2
    kept = []
    for row in rows:
        along = row.period / np.linalg.norm(row.period)
        spans = [
            sorted(
                float(np.dot([point.real, point.imag], along))
                for point in (_placed(turn, int(body[owner]), complex(*end)) for end in pair)
            )
            for owner, pair in zip(row.owners, row.ends, strict=True)
        ]
        low = max(first for first, _ in spans) + reach
        high = min(last for _, last in spans) - reach
        slack = 1e-9 * reach
        kept.append(
            np.array(
                [
                    bool(np.all((s >= low - slack) & (s <= high + slack)))
                    for s in (center[cell] @ along for cell in row.cells)
                ],
                dtype=bool,
            )
        )
    return kept


def solve_rows(
    rows: Sequence[Row],
    mesh: Mesh,
    rest: Mesh,
    body: Index,
    pivots: Array,
    reach: float,
    totals: tuple[Array, Array, Array, Array],
    kept: Sequence[NDArray[np.bool_]] | None = None,
) -> int:
    """Add to ``totals``, capacitance, gradient, Hessian, and element charges per unit
    permittivity, what the cells of ``rows`` hold, equations (154) and (155), and return
    the number of cells solved.

    ``mesh`` is the mesh with the bodies where they are and ``rest`` the same mesh at
    rest; ``body`` gives the body of each conductor and ``pivots`` the point each body
    turns about, as for ``charge_derivatives``. ``kept`` tells for each row which of its
    cells count, by default all. Rows whose cells are translates at rest and lie on the
    same conductors form a family. The counted cell nearest the middle of the first row
    of a family with one, solved with the periodic kernel, gives capacitance and its
    derivatives in travel; a window of ``reach`` around it, with its bodies turned
    rigidly about the cell center, gives those in turn. Every counted cell of the family
    is that cell moved by the turn of its bodies, carried to second order.
    """
    capacitance, gradient, hessian, charge = totals
    turn = _turns_of(mesh, rest, body)
    counted = kept if kept is not None else [np.ones(len(row.cells), dtype=bool) for row in rows]
    families: dict[bytes, list[tuple[Row, NDArray[np.bool_]]]] = {}
    for row, mask in zip(rows, counted, strict=True):
        middle = row.cells[row.middle]
        reference = ((rest.start[middle] + rest.end[middle]) / 2).mean(axis=0)
        order, present, moving, key = canonical(rest, middle, body, reference, {})
        key += present.tobytes() + np.asarray(row.period).tobytes()
        if mask.any():
            families.setdefault(key, []).append((row, mask))
    cache: dict[bytes, tuple[Array, ...]] = {}
    three = np.arange(3)
    for family in families.values():
        first, mask = family[0]
        chosen = np.flatnonzero(mask)
        pick = int(chosen[np.argmin(np.abs(chosen - first.middle))])
        origin = complex(*first.centers[pick])
        spots = {
            b: _placed(turn, b, origin)
            for b in np.unique(body[mesh.conductor[first.cells[pick]]]).tolist()
        }
        found, order, present, moving = _cell(
            mesh, first.cells[pick], body, spots, first, reach, cache
        )
        held, flat, curve, cell_charge = found
        size, rigid = len(present), len(moving)
        flat = flat.reshape(size, size, 3 * rigid)
        curve = curve.reshape(size, size, 3 * rigid, 3 * rigid)
        total = np.zeros((size, size))
        tilt = np.zeros((size, size, 3 * rigid))
        swing = np.zeros((size, size, 3 * rigid, 3 * rigid))
        for row, mask in family:
            for k in np.flatnonzero(mask):
                center = row.centers[k]
                step = np.zeros(3 * rigid)
                arms = np.zeros((rigid, 2))
                offset = complex(*center) - origin
                for j, b in enumerate(moving):
                    shift = (np.exp(1j * turn[int(b)][0]) - 1) * offset
                    step[3 * j : 3 * j + 2] = shift.real, shift.imag
                    arm = _placed(turn, int(b), complex(*center)) - complex(*pivots[b])
                    arms[j] = arm.real, arm.imag
                total += held + flat @ step + 0.5 * (curve @ step) @ step
                moved_slope, moved_bend = repivot(flat + curve @ step, curve, arms)
                tilt += moved_slope
                swing += moved_bend
                charge[np.ix_(row.cells[k][order], present)] = cell_charge
        capacitance[np.ix_(present, present)] += total
        gradient[np.ix_(present, present, moving, three)] += tilt.reshape(size, size, rigid, 3)
        hessian[np.ix_(present, present, moving, three, moving, three)] += swing.reshape(
            size, size, rigid, 3, rigid, 3
        )
    return len(families)


def _turns_of(mesh: Mesh, rest: Mesh, body: Index) -> dict[int, tuple[float, complex, complex]]:
    """Turn of every body from rest, with a point of it at rest and where it is now."""
    owner = body[mesh.conductor]
    turn = {}
    for b in np.unique(owner).tolist():
        on = np.flatnonzero(owner == b)
        k = int(on[np.argmax(np.linalg.norm(rest.end[on] - rest.start[on], axis=1))])
        was = complex(*(rest.end[k] - rest.start[k]))
        now = complex(*(mesh.end[k] - mesh.start[k]))
        turn[b] = (float(np.angle(now / was)), complex(*rest.start[k]), complex(*mesh.start[k]))
    return turn


def _placed(turn: dict[int, tuple[float, complex, complex]], b: int, point: complex) -> complex:
    """Where a point at rest on body ``b`` is now."""
    angle, anchor_rest, anchor_now = turn[b]
    return complex(anchor_now + np.exp(1j * angle) * (point - anchor_rest))


def _cell(
    mesh: Mesh,
    middle: Index,
    body: Index,
    spots: dict[int, complex],
    row: Row,
    reach: float,
    cache: dict[bytes, tuple[Array, ...]],
) -> tuple[tuple[Array, ...], Index, Index, Index]:
    """Capacitance, gradient, Hessian, and element charges of the middle cell in an order
    that congruent cells share, with that order, its conductors, and its bodies; the
    derivatives in turn come from a window around it, equation (155)."""
    from connflex_core.electrostatics.dense import charge_derivatives, solve_charges
    from connflex_core.progress import quiet

    center = (mesh.start + mesh.end) / 2
    reference = center[middle].mean(axis=0)
    order, present, moving, key = canonical(mesh, middle, body, reference, spots)
    key += np.round(np.asarray(row.period) / KEY_GRID).astype(np.int64).tobytes()
    place = {int(g): k for k, g in enumerate(present)}
    ranks = np.array([place[int(g)] for g in mesh.conductor[middle[order]]], dtype=np.intp)
    turning = np.array([[spots[int(b)].real, spots[int(b)].imag] for b in moving])
    kinds = np.array([list(moving).index(b) for b in body[present]])
    if key not in cache:
        cell = Mesh(mesh.start[middle[order]], mesh.end[middle[order]], ranks)
        charges = solve_charges(cell, len(present), quiet, period=row.period)
        slope, bend = charge_derivatives(charges, kinds, turning, quiet)
        cache[key] = (charges.capacitance, slope, bend, charges.charge)
    held, slope, bend, charge = cache[key]
    overlap = max(reach, 2 * float(np.linalg.norm(row.period)))
    low, high = center[middle].min(axis=0) - overlap, center[middle].max(axis=0) + overlap
    inside = np.flatnonzero(np.all((center >= low) & (center <= high), axis=1))
    ordered, near, kinds_near, window_key = canonical(mesh, inside, body, reference, spots)
    window_key += np.isin(inside[ordered], middle).tobytes() + key
    if window_key not in cache:
        elements = inside[ordered]
        within = {int(g): k for k, g in enumerate(near)}
        window = Mesh(
            mesh.start[elements],
            mesh.end[elements],
            np.array([within[int(g)] for g in mesh.conductor[elements]], dtype=np.intp),
        )
        solved = solve_charges(window, len(near), quiet, np.isin(elements, middle))
        bodies_near = [int(b) for b in kinds_near]
        pivot_near = np.array(
            [[spots[b].real, spots[b].imag] if b in spots else list(reference) for b in bodies_near]
        )
        window_slope, window_bend = charge_derivatives(
            solved, np.array([bodies_near.index(int(body[g])) for g in near]), pivot_near, quiet
        )
        mine = np.array([within[int(g)] for g in present])
        theirs = np.array([bodies_near.index(int(b)) for b in moving])
        pick = np.ix_(mine, mine)
        cache[window_key] = (
            window_slope[pick][:, :, theirs],
            window_bend[pick][:, :, theirs][:, :, :, :, theirs],
        )
    window_slope, window_bend = cache[window_key]
    slope, bend = slope.copy(), bend.copy()
    slope[:, :, :, 2] = window_slope[:, :, :, 2]
    bend[:, :, :, 2] = window_bend[:, :, :, 2]
    bend[:, :, :, :, :, 2] = window_bend[:, :, :, :, :, 2]
    return (held, slope, bend, charge), order, present, moving

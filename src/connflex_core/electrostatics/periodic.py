"""Rows of identical cells: the periodic kernel, arrays in a mesh, and their capacitance."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import networkx as nx
import numpy as np
import shapely
from numpy.typing import NDArray
from shapely import MultiPolygon, STRtree

from connflex_core.electrostatics.kernel import pair_derivatives, pair_integrals
from connflex_core.electrostatics.outlines import Mesh

Array = NDArray[np.float64]
Complex = NDArray[np.complex128]
Index = NDArray[np.intp]

LONGEST_MOTIF = 16
GRID = 1e-10
KEY_GRID = 1e-11

IMAGES = 2
_NODES, _WEIGHTS = np.polynomial.legendre.leggauss(4)
_SPOTS = (_NODES + 1) / 2
_SHARES = _WEIGHTS / 2
_SMALL = 0.25


def periodic_integrals(
    a0: Complex, a1: Complex, b0: Complex, b1: Complex, period: complex
) -> Array:
    """Double integral over element pairs of ln|2 sin(pi w / p)| less a constant,
    equation (151), for a row repeated without end along ``period``.

    The nearest ``IMAGES`` copies on either side enter through the exact pair integral,
    the rest through the remainder of (152) by Gauss-Legendre quadrature on both
    elements. A constant per pair is dropped, which a neutral row does not feel.
    """
    a0, a1, b0, b1 = np.broadcast_arrays(a0, a1, b0, b1)
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
    gradient = np.zeros((6, *a0.shape))
    hessian = np.zeros((6, 6, *a0.shape))
    for k in range(-IMAGES, IMAGES + 1):
        slope, bend = pair_derivatives(
            a0,
            a1,
            b0 + k * period,
            b1 + k * period,
            np.asarray(pivot_a, dtype=np.complex128),
            np.asarray(pivot_b + k * period, dtype=np.complex128),
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
        -along * 1j * (second - pivot_b) * np.ones(first.shape),
    ]
    turns = {
        (2, 2): -along * (first - pivot_a) * np.ones(second.shape),
        (5, 5): along * (second - pivot_b) * np.ones(first.shape),
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


def _spots(start: Complex, end: Complex) -> Complex:
    spots: Complex = start[..., None] + (end - start)[..., None] * _SPOTS
    return spots


def _integrate(values: Array) -> Array:
    integral: Array = np.einsum("...ab,a,b->...", values, _SHARES, _SHARES)
    return integral


def _remainder(w: Complex, reach: float) -> tuple[Array, Complex, Complex]:
    """Value, slope, and curvature of the remainder of equation (152) at ``w``, the
    separation of two points turned so that the period lies along the real axis."""
    x = np.pi * w / reach
    small = np.abs(x) < _SMALL
    safe = np.where(small, 1.0, x)
    x2 = np.where(small, x * x, 0.0)
    log_sinc = np.where(
        small,
        -x2 / 6 - x2**2 / 180 - x2**3 / 2835 - x2**4 / 37800,
        np.log(np.sin(safe) / safe),
    )
    cot_less = np.where(
        small,
        -x / 3 - x * x2 / 45 - 2 * x * x2**2 / 945 - x * x2**3 / 4725,
        np.cos(safe) / np.sin(safe) - 1 / safe,
    )
    csc_less = np.where(
        small,
        -1 / 3 - x2 / 15 - 2 * x2**2 / 189 - x2**3 / 675,
        1 / safe**2 - 1 / np.sin(safe) ** 2,
    )
    value: Array = np.asarray(log_sinc.real, dtype=np.float64)
    slope: Complex = np.asarray(np.pi / reach * cot_less, dtype=np.complex128)
    curve: Complex = np.asarray((np.pi / reach) ** 2 * csc_less, dtype=np.complex128)
    for k in range(1, IMAGES + 1):
        far = (k * reach) ** 2
        value = value - np.log(np.abs(1 - w * w / far))
        slope = slope + 2 * w / (far - w * w)
        curve = curve + 2 * (far + w * w) / (far - w * w) ** 2
    return value, slope, curve


@dataclass(frozen=True, eq=False)
class Row:
    """Interior cells of a row of identical cells in a mesh, equation (153).

    Attributes:
        period: Vector from one cell to the next, in meters.
        cells: Elements of each interior cell, in order along ``period``.
        centers: Center of each interior cell at rest, in meters.
        middle: Position in ``cells`` of the cell solved for all of them.
    """

    period: Array
    cells: tuple[Index, ...]
    centers: Array
    middle: int

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
                arrays.append(found_array)
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
                np.lexsort(
                    (
                        (center[members] - k * period) @ across,
                        center[members] @ along - k * pitch,
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


def solve_rows(
    rows: Sequence[Row],
    mesh: Mesh,
    rest: Mesh,
    body: Index,
    pivots: Array,
    reach: float,
    totals: tuple[Array, Array, Array, Array],
) -> int:
    """Add to ``totals``, capacitance, gradient, Hessian, and element charges per unit
    permittivity, what the cells of ``rows`` hold, equations (154) and (155), and return
    the number of cells solved.

    ``mesh`` is the mesh with the bodies where they are and ``rest`` the same mesh at
    rest; ``body`` gives the body of each conductor and ``pivots`` the point each body
    turns about, as for ``charge_derivatives``. Rows whose cells are translates at rest
    and lie on the same conductors form a family. The middle cell of the first row of a
    family, solved with the periodic kernel, gives capacitance and its derivatives in
    travel; a window of ``reach`` around it, with its bodies turned rigidly about the
    cell center, gives those in turn. Every cell of the family is that cell moved by the
    turn of its bodies, carried to second order.
    """
    capacitance, gradient, hessian, charge = totals
    turn = _turns_of(mesh, rest, body)
    families: dict[bytes, list[Row]] = {}
    for row in rows:
        middle = row.cells[row.middle]
        reference = ((rest.start[middle] + rest.end[middle]) / 2).mean(axis=0)
        order, present, moving, key = canonical(rest, middle, body, reference, {})
        key += present.tobytes() + np.asarray(row.period).tobytes()
        families.setdefault(key, []).append(row)
    cache: dict[bytes, tuple[Array, ...]] = {}
    three = np.arange(3)
    for family in families.values():
        first = family[0]
        origin = complex(*first.centers[first.middle])
        spots = {
            b: _placed(turn, b, origin)
            for b in np.unique(body[mesh.conductor[first.cells[first.middle]]]).tolist()
        }
        found, order, present, moving = _cell(
            mesh, first.cells[first.middle], body, spots, first, reach, cache
        )
        held, flat, curve, cell_charge = found
        size, rigid = len(present), len(moving)
        flat = flat.reshape(size, size, 3 * rigid)
        curve = curve.reshape(size, size, 3 * rigid, 3 * rigid)
        total = np.zeros((size, size))
        tilt = np.zeros((size, size, 3 * rigid))
        swing = np.zeros((size, size, 3 * rigid, 3 * rigid))
        for row in family:
            for k, center in enumerate(row.centers):
                step = np.zeros(3 * rigid)
                jacobian = np.eye(3 * rigid)
                lever = np.zeros((rigid, 2))
                offset = complex(*center) - origin
                for j, b in enumerate(moving):
                    shift = (np.exp(1j * turn[int(b)][0]) - 1) * offset
                    step[3 * j : 3 * j + 2] = shift.real, shift.imag
                    arm = _placed(turn, int(b), complex(*center)) - complex(*pivots[b])
                    jacobian[3 * j, 3 * j + 2] = -arm.imag
                    jacobian[3 * j + 1, 3 * j + 2] = arm.real
                    lever[j] = arm.real, arm.imag
                pulled = flat + curve @ step
                total += held + flat @ step + 0.5 * (curve @ step) @ step
                tilt += pulled @ jacobian
                swing += np.einsum("pq,abpr,rs->abqs", jacobian, curve, jacobian)
                for j in range(rigid):
                    swing[:, :, 3 * j + 2, 3 * j + 2] -= (
                        pulled[:, :, 3 * j] * lever[j, 0] + pulled[:, :, 3 * j + 1] * lever[j, 1]
                    )
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


def canonical(
    mesh: Mesh, elements: Index, body: Index, reference: Array, spots: dict[int, complex]
) -> tuple[Index, Index, Index, bytes]:
    """An order of ``elements`` that translates share, their conductors and bodies by
    first appearance in it, and a key of that geometry relative to ``reference``."""
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
        + np.array([list(moving).index(b) for b in owned], dtype=np.int64).tobytes()
        + turned.tobytes()
    )
    return order, present, moving, key

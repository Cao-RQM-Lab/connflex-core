"""Sidewall capacitance from facing walls taken as parallel plates."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from itertools import pairwise
from typing import Any

import numpy as np
import shapely
from numpy.typing import NDArray
from shapely import MultiPolygon, STRtree

from connflex_core.electrostatics.jets import Jet

Array = NDArray[np.float64]
Index = NDArray[np.intp]
Mark = tuple[float, int, int, int]

OWN, PROJECTED, CUTOFF = 0, 1, 2
PIECES_PER_PASS = 8000
SERIES_TERMS = 13
SERIES_BELOW = 0.25
STEEPEST = 10.0
RUN_JUMP = 0.5


@dataclass(frozen=True, eq=False)
class FacingWalls:
    """Pieces of wall that face another conductor, each taken as a plate.

    Attributes:
        start: P by 2 first end points of pieces, in meters.
        end: P by 2 second end points. Going from ``start`` to ``end``, silicon of the
            piece lies to the left and the wall it faces to the right.
        conductor: Conductor of each piece.
        facing: Conductor each piece faces.
        gap: P by 2 distances along the normal of a piece to the wall it faces, at
            ``start`` and at ``end``.
    """

    start: Array
    end: Array
    conductor: Index
    facing: Index
    gap: Array

    @classmethod
    def empty(cls) -> FacingWalls:
        none = np.zeros(0, dtype=np.intp)
        return cls(np.zeros((0, 2)), np.zeros((0, 2)), none, none, np.zeros((0, 2)))


@dataclass(frozen=True, eq=False)
class FacingRuns:
    """Pieces of wall that face another conductor, with their length and its
    derivatives for rigid motion of the bodies that set it, equation (103).

    Attributes:
        edge: Outline edge of each piece, numbered as by ``outline_edges``.
        conductor: Conductor of each piece.
        facing: Conductor each piece faces.
        start: P by 2 first end points, as on ``FacingWalls``.
        end: P by 2 second end points.
        gap: P by 2 gaps at ``start`` and at ``end``.
        bodies: P by 4 bodies of the parameters of ``length``: of the edge, of the edge
            it faces, and of the edges that set its two ends, -1 where fewer.
        length: Length of each piece, with its gradient and Hessian in travel along x
            and y and turn about its pivot of each of the four bodies.
        distance: Gap at the middle of each piece along the normal of its wall, with
            its derivatives in the same parameters.
        middle: Distance of that middle along its edge from the edge's first point,
            with its derivatives.
    """

    edge: Index
    conductor: Index
    facing: Index
    start: Array
    end: Array
    gap: Array
    bodies: Index
    length: Jet
    distance: Jet
    middle: Jet


def facing_plates(
    regions: Sequence[MultiPolygon], body: Index, pivots: Array, reach: float
) -> tuple[Array, Array, Array, FacingWalls]:
    """Capacitance per permittivity between facing walls, equations (97) to (103).

    ``body`` gives the rigid body of each region and ``pivots`` the point each body
    turns about. A wall counts toward two conductors when the normal from one first
    meets the other within ``reach``. Returns the n by n matrix with rows that sum to
    zero, its gradient and Hessian in the layout of ``SidewallCapacitance``, and the
    pieces of wall.
    """
    edges, pieces, slot, used, members = _pieces(regions, body, reach)
    start, end, conductor, outward, length = edges
    moving = body[conductor]
    count, rigid = len(regions), len(pivots)
    between = np.zeros((count, count))
    slopes = np.zeros((count, count, 3 * rigid))
    bends = np.zeros((count, count, 3 * rigid, 3 * rigid))
    left, right = np.zeros(len(pieces)), np.zeros(len(pieces))
    gap = np.zeros((len(pieces), 2))
    three = np.arange(3)
    for slots in range(1, 5):
        group = np.flatnonzero(used == slots)
        for first in range(0, len(group), PIECES_PER_PASS):
            part = group[first : first + PIECES_PER_PASS]
            value, low, high, _, gap[part] = _integrals(
                pieces[part], start, end, length, outward, pivots[moving], slot[part], slots, reach
            )
            left[part], right[part] = low.value, high.value
            near, far = conductor[pieces[part, 0]], conductor[pieces[part, 1]]
            for row, column in ((near, far), (far, near)):
                np.add.at(between, (row, column), value.value / 2)
                if slots == 1:
                    continue
                for i in range(slots):
                    one = 3 * members[part, i, None] + three
                    np.add.at(
                        slopes,
                        (row[:, None], column[:, None], one),
                        value.gradient[:, 3 * i : 3 * i + 3] / 2,
                    )
                    for j in range(slots):
                        other = 3 * members[part, j, None] + three
                        np.add.at(
                            bends,
                            (
                                row[:, None, None],
                                column[:, None, None],
                                one[:, :, None],
                                other[:, None, :],
                            ),
                            value.hessian[:, 3 * i : 3 * i + 3, 3 * j : 3 * j + 3] / 2,
                        )
    walls = _walls(edges, pieces, left, right, gap)
    return (
        maxwell(between),
        maxwell(slopes).reshape(count, count, rigid, 3),
        maxwell(bends).reshape(count, count, rigid, 3, rigid, 3),
        walls,
    )


def facing_runs(
    regions: Sequence[MultiPolygon],
    body: Index,
    pivots: Array,
    reach: float,
    smooth: float = 0.0,
) -> FacingRuns:
    """Pieces of wall that face another conductor within ``reach``, as for
    ``facing_plates``, with the length of each and its derivatives.

    With ``smooth``, a piece whose ends are the ends of its edge or of the edge it faces
    takes the overlap of the two rounded over that length where ends pass each other,
    equation (161), so that its length has continuous slope as walls slide.
    """
    edges, pieces, slot, used, members = _pieces(regions, body, reach)
    start, end, conductor, outward, length = edges
    moving = body[conductor]
    count = len(pieces)
    left, right = np.zeros(count), np.zeros(count)
    gap = np.zeros((count, 2))
    run, distance, halfway = (
        Jet(np.zeros(count), np.zeros((count, 12)), np.zeros((count, 12, 12))) for _ in range(3)
    )
    for slots in range(1, 5):
        group = np.flatnonzero(used == slots)
        for first in range(0, len(group), PIECES_PER_PASS):
            part = group[first : first + PIECES_PER_PASS]
            _, low, high, middle, gap[part], (u0, u1) = _integrals(
                pieces[part],
                start,
                end,
                length,
                outward,
                pivots[moving],
                slot[part],
                slots,
                reach,
                ends=True,
            )
            left[part], right[part] = low.value, high.value
            span = high - low
            if smooth > 0:
                span = _pick(
                    _geometric(pieces[part]),
                    _overlap(u0, u1, length[pieces[part, 0]], smooth),
                    span,
                )
            for into, jet in ((run, span), (distance, middle), (halfway, (low + high) * 0.5)):
                into.value[part] = jet.value
                into.gradient[part, : 3 * slots] = jet.gradient
                into.hessian[part, : 3 * slots, : 3 * slots] = jet.hessian
    walls = _walls(edges, pieces, left, right, gap)
    bodies = np.where(np.arange(4) < used[:, None], members, -1)
    return FacingRuns(
        edge=pieces[:, 0],
        conductor=walls.conductor,
        facing=walls.facing,
        start=walls.start,
        end=walls.end,
        gap=walls.gap,
        bodies=bodies,
        length=run,
        distance=distance,
        middle=halfway,
    )


def _pieces(
    regions: Sequence[MultiPolygon], body: Index, reach: float
) -> tuple[tuple[Array, Array, Index, Array, Array], Index, Index, Index, Index]:
    """Outline edges with their lengths, pieces of them that face another conductor, and
    for each piece the slot of each of its four bodies, how many slots it uses, and the
    body in each slot."""
    start, end, conductor, outward = outline_edges(regions)
    along = end - start
    length = np.hypot(along[:, 0], along[:, 1])
    tangent = along / length[:, None]
    normal = outward[:, None] * np.stack([tangent[:, 1], -tangent[:, 0]], axis=1)
    pieces = _visible(start, end, conductor, length, tangent, normal, reach)

    # Bodies of a piece: of its edge, of the edge it faces, and of the edges that set
    # its two ends. Each distinct one takes a slot of three parameters.
    moving = body[conductor]
    involved = moving[pieces[:, [0, 1, 3, 6]]]
    each = np.arange(len(pieces))
    slot = np.zeros((len(pieces), 4), dtype=np.intp)
    used = np.ones(len(pieces), dtype=np.intp)
    members = np.zeros((len(pieces), 4), dtype=np.intp)
    members[:, 0] = involved[:, 0]
    for k in range(1, 4):
        match = involved[:, :k] == involved[:, k : k + 1]
        known = match.any(axis=1)
        slot[:, k] = np.where(known, slot[each, match.argmax(axis=1)], used)
        members[each, slot[:, k]] = involved[:, k]
        used += ~known
    return (start, end, conductor, outward, length), pieces, slot, used, members


def _walls(
    edges: tuple[Array, Array, Index, Array, Array],
    pieces: Index,
    left: Array,
    right: Array,
    gap: Array,
) -> FacingWalls:
    """Pieces as ``FacingWalls`` from their ends along their edges."""
    start, end, conductor, outward, length = edges
    tangent = (end - start) / length[:, None]
    edge = pieces[:, 0]
    ends = np.stack([left, right], axis=1)
    gap = gap.copy()
    turned = outward[edge] < 0
    ends[turned], gap[turned] = ends[turned, ::-1], gap[turned, ::-1]
    return FacingWalls(
        start=start[edge] + ends[:, :1] * tangent[edge],
        end=start[edge] + ends[:, 1:] * tangent[edge],
        conductor=conductor[edge],
        facing=conductor[pieces[:, 1]],
        gap=gap,
    )


def maxwell(between: Array) -> Array:
    """Capacitance matrix from capacitances between pairs, equation (108): negated, with
    the row sums on the diagonal."""
    matrix = -between
    diagonal = np.arange(len(between))
    matrix[diagonal, diagonal] = between.sum(axis=1)
    return matrix


def outline_edges(regions: Sequence[MultiPolygon]) -> tuple[Array, Array, Index, Array]:
    """Start, end, conductor, and side of the outward normal of every outline edge."""
    starts, ends, owners, sides = [], [], [], []
    for index, region in enumerate(regions):
        rings, polygon_of = shapely.get_rings(shapely.get_parts(region), return_index=True)
        coords, ring_of = shapely.get_coordinates(rings, return_index=True)
        joined = ring_of[:-1] == ring_of[1:]
        a, b, ring = coords[:-1][joined], coords[1:][joined], ring_of[:-1][joined]
        exterior = np.ones(rings.size, dtype=bool)
        exterior[1:] = polygon_of[1:] != polygon_of[:-1]
        twice_area = np.bincount(
            ring, weights=a[:, 0] * b[:, 1] - b[:, 0] * a[:, 1], minlength=rings.size
        )
        side = np.where(exterior == (twice_area > 0), 1.0, -1.0)[ring]
        kept = np.any(a != b, axis=1)
        starts.append(a[kept])
        ends.append(b[kept])
        sides.append(side[kept])
        owners.append(np.full(int(kept.sum()), index, dtype=np.intp))
    return (
        np.concatenate(starts),
        np.concatenate(ends),
        np.concatenate(owners),
        np.concatenate(sides),
    )


def _visible(
    start: Array,
    end: Array,
    conductor: Index,
    length: Array,
    tangent: Array,
    normal: Array,
    reach: float,
) -> Index:
    """Pieces of edges that face another conductor, equations (99) to (101).

    Each row holds the edge, the edge it faces, and for each end of the piece its kind,
    the edge that sets it, and which end point of that edge. An edge that rises more
    than ``STEEPEST`` times its run along a wall hides what lies behind it and is no
    plate with the wall.
    """
    lines = np.asarray(shapely.linestrings(np.stack([start, end], axis=1)))
    swept = shapely.polygons(
        np.stack([start, end, end + reach * normal, start + reach * normal], axis=1)
    )
    first, second = STRtree(lines).query(swept, predicate="intersects")
    apart = first != second
    first, second = first[apart], second[apart]
    u = np.stack(
        [np.sum((p[second] - start[first]) * tangent[first], axis=1) for p in (start, end)]
    )
    v = np.stack([np.sum((p[second] - start[first]) * normal[first], axis=1) for p in (start, end)])
    tiny = 1e-9 * length[first]
    span = u[1] - u[0]
    crossing = np.abs(span) > tiny
    slope = np.where(crossing, (v[1] - v[0]) / np.where(crossing, span, 1.0), 0.0)
    low, high = np.maximum(u.min(axis=0), 0.0), np.minimum(u.max(axis=0), length[first])
    height = v[0] + slope * ((low + high) / 2 - u[0])
    seen = crossing & (high - low > tiny) & (height > 1e-9 * reach)
    level = np.abs(v[1] - v[0]) <= STEEPEST * np.abs(span)
    order = np.lexsort((low[seen], first[seen]))
    columns = (first, second, low, high, u[0], v[0], slope, u.min(axis=0), u.max(axis=0), level)
    edges, others, lows, highs, u0s, v0s, slopes, firsts, lasts, levels = (
        column[seen][order].tolist() for column in columns
    )
    lower_end = np.argmin(u, axis=0)[seen][order].tolist()
    lengths, owners = length.tolist(), conductor.tolist()

    pieces: list[list[int]] = []
    position, total = 0, len(edges)
    while position < total:
        edge = edges[position]
        stop = position
        while stop < total and edges[stop] == edge:
            stop += 1
        full = lengths[edge]
        tolerance = 1e-9 * full
        marks: list[Mark] = [(0.0, OWN, edge, 0), (full, OWN, edge, 1)]
        for k in range(position, stop):
            if firsts[k] > tolerance:
                marks.append((lows[k], PROJECTED, others[k], lower_end[k]))
            if lasts[k] < full - tolerance:
                marks.append((highs[k], PROJECTED, others[k], 1 - lower_end[k]))
            if slopes[k] != 0.0:
                cut = u0s[k] + (reach - v0s[k]) / slopes[k]
                if lows[k] < cut < highs[k]:
                    marks.append((cut, CUTOFF, others[k], 0))
        marks.sort()
        groups = [[marks[0]]]
        for mark in marks[1:]:
            if mark[0] - groups[-1][0][0] > tolerance:
                groups.append([mark])
            else:
                groups[-1].append(mark)

        nearest: list[int] = []
        faced: list[int] = []
        active: list[int] = []
        upcoming = position
        for before, after in pairwise(groups):
            middle = (before[0][0] + after[0][0]) / 2
            while upcoming < stop and lows[upcoming] <= middle:
                active.append(upcoming)
                upcoming += 1
            active = [k for k in active if highs[k] >= middle]
            best, lowest, plate = -1, reach, False
            for k in active:
                distance = v0s[k] + slopes[k] * (middle - u0s[k])
                if distance < lowest:
                    best, lowest, plate = others[k], distance, levels[k]
            nearest.append(best)
            faced.append(best if plate and owners[best] != owners[edge] else -1)

        k = 0
        while k < len(faced):
            best, last = faced[k], k
            if best >= 0:
                while last + 1 < len(faced) and faced[last + 1] == best:
                    last += 1
                previous = nearest[k - 1] if k else -1
                following = nearest[last + 1] if last + 1 < len(faced) else -1
                opening = _setting(groups[k], best, previous)
                closing = _setting(groups[last + 1], best, following)
                pieces.append([edge, best, *opening[1:], *closing[1:]])
            k = last + 1
        position = stop
    return np.array(pieces, dtype=np.intp).reshape(-1, 8)


def _setting(marks: list[Mark], facing: int, neighbor: int) -> Mark:
    """Mark that sets an end of a piece where several coincide.

    An end of the facing edge comes first, then an end of the edge itself, then the
    edge faced on the other side of the mark.
    """

    def rank(mark: Mark) -> int:
        if mark[2] == facing:
            return 0
        if mark[1] == OWN:
            return 1
        return 2 if mark[2] == neighbor else 3

    return min(marks, key=rank)


def _integrals(
    pieces: Index,
    start: Array,
    end: Array,
    length: Array,
    outward: Array,
    pivot: Array,
    slot: Index,
    slots: int,
    reach: float,
    ends: bool = False,
) -> tuple[Any, ...]:
    """Integral of inverse gap over each piece with its derivatives, equations (102) and (103).

    Parameters are travel along x and y and rotation about its pivot for each of
    ``slots`` bodies. Columns of ``slot`` give the body of the edge, of the edge it
    faces, and of the edges that set the two ends of the piece; ``pivot`` holds the
    pivot of every edge. Also returns both ends of each piece as distances along its
    edge and the gap at its middle, with their derivatives, and the gap at both ends;
    with ``ends`` as well where both ends of the faced edge fall along the edge.
    """
    edge, other = pieces[:, 0], pieces[:, 1]

    def moved(points: Array, of: Index, column: int) -> tuple[Jet, Jet]:
        return _moved(points, pivot[of], slot[:, column], 3 * slots)

    p0, p1 = moved(start[edge], edge, 0), moved(end[edge], edge, 0)
    span = length[edge]
    tangent = ((p1[0] - p0[0]) / span, (p1[1] - p0[1]) / span)
    normal = (tangent[1] * outward[edge], tangent[0] * -outward[edge])

    def frame(points: Array, of: Index, column: int) -> tuple[Jet, Jet]:
        x, y = moved(points, of, column)
        dx, dy = x - p0[0], y - p0[1]
        return dx * tangent[0] + dy * tangent[1], dx * normal[0] + dy * normal[1]

    u0, v0 = frame(start[other], other, 1)
    u1, v1 = frame(end[other], other, 1)
    slope = (v1 - v0) / (u1 - u0)

    def mark(column: int) -> Jet:
        kind, setter, which = (
            np.asarray(pieces[:, 3 * column - 4 + k], dtype=np.intp) for k in range(3)
        )
        point = np.where(which[:, None] == 0, start[setter], end[setter])
        projected = frame(point, setter, column)[0]
        own, cast, sloped = (
            np.asarray(kind == one, dtype=np.float64) for one in (OWN, PROJECTED, CUTOFF)
        )
        cut = u0 + (reach - v0) / (slope * sloped + (1 - sloped))
        return Jet.constant(which * span * own, 3 * slots) + projected * cast + cut * sloped

    left, right = mark(2), mark(3)
    run = right - left
    gap = v0 + slope * ((left + right) * 0.5 - u0)
    value = run / gap * (slope * run / (gap * 2.0)).through(_shape, _shape_slope, _shape_bend)
    gaps = np.stack([(v0 + slope * (at - u0)).value for at in (left, right)], axis=1)
    if ends:
        return value, left, right, gap, gaps, (u0, u1)
    return value, left, right, gap, gaps


def _geometric(pieces: Index) -> NDArray[np.bool_]:
    """Pieces whose two ends are ends of their own edge or of the edge they face."""
    edge, other = pieces[:, 0], pieces[:, 1]
    found = np.ones(len(pieces), dtype=bool)
    for kind, setter in ((pieces[:, 2], pieces[:, 3]), (pieces[:, 5], pieces[:, 6])):
        found &= ((kind == OWN) & (setter == edge)) | ((kind == PROJECTED) & (setter == other))
    return found


def _overlap(first: Jet, second: Jet, span: Array, width: float) -> Jet:
    """Overlap of an edge from 0 to ``span`` with the faced edge from ``first`` to
    ``second`` along it, rounded over ``width`` where ends pass, equation (161)."""
    low = _pick(first.value <= second.value, first, second)
    high = _pick(first.value <= second.value, second, first)

    def onset(x: Jet) -> Jet:
        return x.through(
            lambda u: _onset(u, width),
            lambda u: _onset_slope(u, width),
            lambda u: _onset_bend(u, width),
        )

    def even(x: Jet) -> Jet:
        return x.through(
            lambda u: _even(u, width),
            lambda u: _even_slope(u, width),
            lambda u: _even_bend(u, width),
        )

    # Each of the four ramps of the overlap is rounded on the side where the walls still
    # face each other, the two where they stop facing so that the overlap meets zero
    # with zero slope.
    return onset(span - low) - even(span - high) - even(-low) - high + onset(high)


def _onset(x: Array, width: float) -> Array:
    """Zero below 0, x above ``width``, and between them a cubic with matching slopes:
    where walls begin to face each other."""
    inner = np.clip(x, 0.0, width)
    return np.where(x >= width, x, 2 * inner**2 / width - inner**3 / width**2)


def _onset_slope(x: Array, width: float) -> Array:
    inner = np.clip(x, 0.0, width)
    return np.where(x >= width, 1.0, 4 * inner / width - 3 * inner**2 / width**2)


def _onset_bend(x: Array, width: float) -> Array:
    inner = np.clip(x, 0.0, width)
    return np.where((x > 0) & (x < width), 4 / width - 6 * inner / width**2, 0.0)


def _even(x: Array, width: float) -> Array:
    """Zero below and x above half ``width`` either side of 0, and a parabola between:
    where an end of one wall passes an end of the other."""
    inner = np.clip(x, -width / 2, width / 2)
    return np.where(x >= width / 2, x, (inner + width / 2) ** 2 / (2 * width))


def _even_slope(x: Array, width: float) -> Array:
    inner = np.clip(x, -width / 2, width / 2)
    return np.where(x >= width / 2, 1.0, (inner + width / 2) / width)


def _even_bend(x: Array, width: float) -> Array:
    return np.where(abs(x) < width / 2, 1 / width, 0.0)


def _pick(mask: NDArray[np.bool_], yes: Jet, no: Jet) -> Jet:
    """Jet of ``yes`` where ``mask`` holds and of ``no`` elsewhere."""
    return Jet(
        np.where(mask, yes.value, no.value),
        np.where(mask[:, None], yes.gradient, no.gradient),
        np.where(mask[:, None, None], yes.hessian, no.hessian),
    )


def _moved(points: Array, pivot: Array, slot: Index, parameters: int) -> tuple[Jet, Jet]:
    """Coordinates of points that each move rigidly with one body, equation (83)."""
    arm = points - pivot
    rows = np.arange(len(points))
    result = []
    for axis, turned in ((0, -arm[:, 1]), (1, arm[:, 0])):
        gradient = np.zeros((len(points), parameters))
        gradient[rows, 3 * slot + axis] = 1.0
        gradient[rows, 3 * slot + 2] = turned
        hessian = np.zeros((len(points), parameters, parameters))
        hessian[rows, 3 * slot + 2, 3 * slot + 2] = -arm[:, axis]
        result.append(Jet(points[:, axis].copy(), gradient, hessian))
    return result[0], result[1]


def _shape(x: Array) -> Array:
    """atanh(x) / x, the second factor of equation (102) for walls that are not parallel."""
    small = np.abs(x) < SERIES_BELOW
    safe = np.where(small, 0.5, x)
    series = sum(x ** (2 * k) / (2 * k + 1) for k in range(SERIES_TERMS))
    return np.where(small, series, np.arctanh(safe) / safe)


def _shape_slope(x: Array) -> Array:
    small = np.abs(x) < SERIES_BELOW
    safe = np.where(small, 0.5, x)
    series = sum(2 * k * x ** (2 * k - 1) / (2 * k + 1) for k in range(1, SERIES_TERMS))
    return np.where(small, series, (1 / (1 - safe**2) - np.arctanh(safe) / safe) / safe)


def _shape_bend(x: Array) -> Array:
    small = np.abs(x) < SERIES_BELOW
    safe = np.where(small, 0.5, x)
    series = sum(
        2 * k * (2 * k - 1) * x ** (2 * k - 2) / (2 * k + 1) for k in range(1, SERIES_TERMS)
    )
    first = (1 / (1 - safe**2) - np.arctanh(safe) / safe) / safe
    return np.where(small, series, (2 * safe / (1 - safe**2) ** 2 - 2 * first) / safe)


def feature_sizes(walls: FacingWalls, net: Index) -> dict[tuple[int, int], float]:
    """Feature size of each pair of nets that face each other, equation (105).

    ``net`` gives the net of every conductor. Pieces that continue each other along one
    outline, face one net, and meet it where the last piece left off form one run. A
    run's length over its mean gap is its capacitance per permittivity, and the feature
    size is the mean of that over the runs between two nets, weighted by itself. Keys
    are pairs of nets, smaller first.
    """
    count = len(walls.conductor)
    if not count:
        return {}
    along = walls.end - walls.start
    length = np.hypot(along[:, 0], along[:, 1])
    right = np.stack([along[:, 1], -along[:, 0]], axis=1) / length[:, None]
    weight = length * _mean_inverse(walls.gap)
    own, faced = net[walls.conductor], net[walls.facing]
    hit = [walls.start + right * walls.gap[:, :1], walls.end + right * walls.gap[:, 1:]]
    grid = 1e-9 * float(np.max(length))
    starts: dict[tuple[int, int, int, int], list[int]] = {}
    for k, (x, y) in enumerate(np.round(walls.start / grid).astype(np.int64).tolist()):
        starts.setdefault((x, y, int(own[k]), int(faced[k])), []).append(k)
    parent = list(range(count))

    def root(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    ends = np.round(walls.end / grid).astype(np.int64).tolist()
    for k, (x, y) in enumerate(ends):
        for j in starts.get((x, y, int(own[k]), int(faced[k])), []):
            jump = float(np.hypot(*(hit[1][k] - hit[0][j])))
            if j != k and jump <= RUN_JUMP * min(walls.gap[k, 1], walls.gap[j, 0]):
                parent[root(j)] = root(k)
    runs: dict[int, list[float]] = {}
    for k in range(count):
        run = runs.setdefault(root(k), [0.0, float(own[k]), float(faced[k])])
        run[0] += weight[k]
    totals: dict[tuple[int, int], list[float]] = {}
    for run_weight, a, b in runs.values():
        pair = (int(min(a, b)), int(max(a, b)))
        total = totals.setdefault(pair, [0.0, 0.0])
        total[0] += run_weight**2
        total[1] += run_weight
    return {pair: aspect / weight for pair, (aspect, weight) in totals.items()}


def _mean_inverse(gap: Array) -> Array:
    """Mean of 1 / gap over a piece whose gap changes linearly between its two ends."""
    low, high = gap.min(axis=1), gap.max(axis=1)
    close = high - low <= 1e-9 * high
    spread = np.where(close, 1.0, high - low)
    mean: Array = np.where(close, 1 / high, np.log(high / np.where(close, 1.0, low)) / spread)
    return mean

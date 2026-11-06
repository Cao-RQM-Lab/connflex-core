"""Fringe of the field above and below the device layer, equations (157) to (163)."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field, replace

import numpy as np
import shapely
from numpy.typing import NDArray
from shapely import MultiPolygon, Polygon, STRtree

from connflex_core.electrostatics.facing import facing_runs, outline_edges
from connflex_core.electrostatics.jets import Jet
from connflex_core.electrostatics.section import Cut, Section, solve_section
from connflex_core.electrostatics.translates import KEY_GRID

Array = NDArray[np.float64]
Index = NDArray[np.intp]

KEY_STEPS = 16
ITEMS_PER_PASS = 400


@dataclass(frozen=True, eq=False)
class FringePlan:
    """What the fringe of a group of device electrodes needs beyond their outlines.

    Attributes:
        handle: Regions of the handle electrodes that may lie beneath.
        handle_body: Body of each handle electrode, numbered among the bodies of the
            group, after those of its device electrodes when it moves with none of them.
        thickness: Thickness of the device layer, in meters.
        box: Height of the released gap above the handle, in meters.
        reach: How far a cut reaches either side of its wall, in meters.
        grading: Shortening of cross-section elements toward corners.
        sections: Solutions met so far, by key of their rounded cut.
        last: The motion last evaluated and its result.
    """

    handle: tuple[MultiPolygon, ...]
    handle_body: Index
    thickness: float
    box: float
    reach: float
    grading: float
    sections: dict[bytes, Section] = field(default_factory=dict)
    last: list[tuple[Array, FringeCapacitance]] = field(default_factory=list)

    @property
    def step(self) -> float:
        """Grid that cuts are rounded to for their key, equation (160)."""
        return min(self.thickness, self.box) / KEY_STEPS


@dataclass(frozen=True, eq=False)
class FringeCapacitance:
    """Capacitance from the fringe above and below the device layer, per unit
    permittivity, equation (163).

    Conductors are the device electrodes of the group, then the handle electrodes.
    Each body has six parameters: travel along x and y and turn about z at device
    mid-plane about its pivot, as for ``SidewallCapacitance``, then rise and turns
    about x and y as in equation (12).

    Attributes:
        capacitance: k by k matrix with rows that sum to zero.
        gradient: k by k by m by 6 derivatives.
        hessian: k by k by m by 6 by m by 6 second derivatives.
        sections: Number of distinct cross-sections the cuts took.
    """

    capacitance: Array
    gradient: Array
    hessian: Array
    sections: int = 0


@dataclass(frozen=True, eq=False)
class _Cutting:
    """One cut: its key and rounded geometry, the conductor and body of the group behind
    each of its labels, and for each of its bodies where its wall nearest the cut wall
    lies and how far that is from the rounding."""

    key: bytes
    cut: Cut
    conductors: Index
    bodies: Index
    position: Array
    offset: Array
    faced: int


@dataclass(frozen=True, eq=False)
class _Items:
    """Cuts with the weight each counts with: the core conductor, where the cut is, its
    normal, the weight, the bodies of the weight's parameters, the gap to the faced
    wall when the exact one is known, and where the cut stands along its edge at rest:
    the edge's first point and direction there and the distance along it."""

    cuttings: list[_Cutting]
    owner: Index
    point: Array
    normal: Array
    weight: Jet
    slots: Index
    distance: Jet | None
    start: Array
    tangent: Array
    middle: Jet


def find_fringe(
    plan: FringePlan,
    regions: Sequence[MultiPolygon],
    body: Index,
    pivots: Array,
    centers: Array,
    motion: Array,
) -> FringeCapacitance:
    """Fringe capacitance of device electrodes ``regions`` where they are, equations
    (157) to (163).

    ``body`` gives the body of each region. ``pivots`` holds where each body of the
    group turns about in the plane now, ``centers`` its center at rest, and ``motion``
    its six parameters from rest, as for ``FringeCapacitance``.
    """
    if plan.last and np.array_equal(plan.last[0][0], motion):
        return plan.last[0][1]
    devices, bodies = len(regions), len(pivots)
    count = devices + len(plan.handle)
    runs = facing_runs(regions, body, pivots, plan.reach, plan.box)
    start, end, conductor, outward = outline_edges(regions)
    along = end - start
    length = np.hypot(along[:, 0], along[:, 1])
    tangent = along / length[:, None]
    normal = outward[:, None] * np.stack([tangent[:, 1], -tangent[:, 0]], axis=1)
    middle = (start + end) / 2

    # Edges carried back to rest with the body of their conductor.
    owner_body = body[conductor]
    back = -motion[owner_body, 2]
    turn = np.stack([[np.cos(back), -np.sin(back)], [np.sin(back), np.cos(back)]]).transpose(
        2, 0, 1
    )
    rest_start = centers[owner_body] + np.einsum("eij,ej->ei", turn, start - pivots[owner_body])
    rest_tangent = np.einsum("eij,ej->ei", turn, tangent)
    point = (runs.start + runs.end) / 2
    facing_normal = normal[runs.edge]
    across = point + facing_normal * runs.gap.mean(axis=1)[:, None] / 2
    owners, polygons = _solid(regions)
    piece_cuts = _cuts(
        plan,
        polygons,
        owners,
        body,
        runs.conductor,
        runs.facing,
        point,
        facing_normal,
        _beneath(plan, across),
        devices,
    )
    edge_cuts = _cuts(
        plan,
        polygons,
        owners,
        body,
        conductor,
        None,
        middle,
        normal,
        _beneath(plan, middle),
        devices,
    )
    for cutting in [*piece_cuts, *edge_cuts]:
        if cutting is not None and cutting.key not in plan.sections:
            plan.sections[cutting.key] = solve_section(
                cutting.cut, plan.thickness, plan.box, plan.reach, plan.grading
            )

    taken = np.flatnonzero([cutting is not None for cutting in piece_cuts])
    open_edge = np.array([cutting is not None for cutting in edge_cuts], dtype=bool)
    covered = taken[open_edge[runs.edge[taken]]]
    lone = np.flatnonzero(open_edge)
    batches = [
        _Items(
            _present(piece_cuts, taken),
            runs.conductor[taken],
            point[taken],
            facing_normal[taken],
            _rows(runs.length, taken),
            runs.bodies[taken],
            _rows(runs.distance, taken),
            rest_start[runs.edge[taken]],
            rest_tangent[runs.edge[taken]],
            _rows(runs.middle, taken),
        ),
        _Items(
            _present(edge_cuts, runs.edge[covered]),
            runs.conductor[covered],
            middle[runs.edge[covered]],
            facing_normal[covered],
            -_rows(runs.length, covered),
            runs.bodies[covered],
            None,
            rest_start[runs.edge[covered]],
            rest_tangent[runs.edge[covered]],
            _rows(runs.middle, covered),
        ),
        _Items(
            _present(edge_cuts, lone),
            conductor[lone],
            middle[lone],
            normal[lone],
            Jet.constant(length[lone], 12),
            np.full((len(lone), 4), -1, dtype=np.intp),
            None,
            rest_start[lone],
            rest_tangent[lone],
            Jet.constant(length[lone] / 2, 12),
        ),
    ]
    total = np.zeros((count, count))
    slope = np.zeros((count, count, 6 * bodies + 1))
    bend = np.zeros((count, count, 6 * bodies + 1, 6 * bodies + 1))
    for batch in batches:
        for low in range(0, len(batch.cuttings), ITEMS_PER_PASS):
            chosen = np.arange(low, min(low + ITEMS_PER_PASS, len(batch.cuttings)))
            part = _Items(
                [batch.cuttings[k] for k in chosen],
                batch.owner[chosen],
                batch.point[chosen],
                batch.normal[chosen],
                _rows(batch.weight, chosen),
                batch.slots[chosen],
                None if batch.distance is None else _rows(batch.distance, chosen),
                batch.start[chosen],
                batch.tangent[chosen],
                _rows(batch.middle, chosen),
            )
            _assemble(plan, part, pivots, centers, motion, total, slope, bend)
    found = _symmetric(total, slope[:, :, :-1], bend[:, :, :-1, :-1], devices, bodies)
    keys = {cutting.key for cutting in [*piece_cuts, *edge_cuts] if cutting is not None}
    found = replace(found, sections=len(keys))
    plan.last[:] = [(motion.copy(), found)]
    return found


def _present(cuttings: list[_Cutting | None], chosen: Index) -> list[_Cutting]:
    """Cuts at ``chosen``, every one of them taken."""
    found = [cuttings[k] for k in chosen.tolist()]
    return [cutting for cutting in found if cutting is not None]


def _rows(jet: Jet, rows: Index) -> Jet:
    return Jet(jet.value[rows], jet.gradient[rows], jet.hessian[rows])


def _beneath(plan: FringePlan, points: Array) -> Index:
    """Handle electrode under each point, -1 where none."""
    found = np.full(len(points), -1, dtype=np.intp)
    if not plan.handle or not len(points):
        return found
    where, which = STRtree(list(plan.handle)).query(shapely.points(points), predicate="intersects")
    found[where] = which
    return found


def _solid(regions: Sequence[MultiPolygon]) -> tuple[Index, list[Polygon]]:
    """Polygons of the regions, with holes that hold no other polygon filled as the
    mesh leaves them, and the region of each."""
    polygons = [polygon for region in regions for polygon in region.geoms]
    owners = np.array([k for k, region in enumerate(regions) for _ in region.geoms], dtype=np.intp)
    tree = STRtree(polygons)
    solid = []
    for position, polygon in enumerate(polygons):
        kept = [
            ring
            for ring in polygon.interiors
            if np.any(tree.query(Polygon(ring), predicate="intersects") != position)
        ]
        solid.append(Polygon(polygon.exterior, kept))
    return owners, solid


def _cuts(
    plan: FringePlan,
    polygons: list[Polygon],
    owners: Index,
    body: Index,
    wall: Index,
    faced: Index | None,
    points: Array,
    normals: Array,
    beneath: Index,
    devices: int,
) -> list[_Cutting | None]:
    """Cut across the wall of conductor ``wall`` at each of ``points``, equation (157).

    A cut reaches ``plan.reach`` behind the wall and, for a piece that faces conductor
    ``faced``, as far ahead of it. One that does not reach ahead holds the wall's own
    side alone and is taken only where the handle lies beneath. None where no cut is
    taken.
    """
    reach = plan.reach
    ahead = faced is not None
    found: list[_Cutting | None] = [None] * len(points)
    wanted = np.arange(len(points)) if ahead else np.flatnonzero(beneath >= 0)
    if not len(wanted):
        return found
    origin, direction = points[wanted], normals[wanted]
    tips = origin + direction * (reach if ahead else 0.0)
    lines = np.asarray(shapely.linestrings(np.stack([origin - direction * reach, tips], axis=1)))
    line_of, polygon_of = STRtree(polygons).query(lines, predicate="intersects")
    crossed = shapely.intersection(lines[line_of], np.asarray(polygons, dtype=object)[polygon_of])
    parts, part_of = shapely.get_parts(crossed, return_index=True)
    coords, coord_of = shapely.get_coordinates(parts, return_index=True)
    line = line_of[part_of][coord_of]
    along = np.einsum("ij,ij->i", coords - origin[line], direction[line])
    low = np.full(len(parts), np.inf)
    high = np.full(len(parts), -np.inf)
    np.minimum.at(low, coord_of, along)
    np.maximum.at(high, coord_of, along)
    tolerance = 1e-6 * reach
    keep = np.flatnonzero(high - low > tolerance)
    per_line: dict[int, list[list[float]]] = {}
    for a, b, which, owner in zip(
        low[keep].tolist(),
        high[keep].tolist(),
        line_of[part_of[keep]].tolist(),
        owners[polygon_of[part_of[keep]]].tolist(),
        strict=True,
    ):
        per_line.setdefault(which, []).append([a, b, owner])
    for position, spans in per_line.items():
        index = int(wanted[position])
        own = int(wall[index])
        spans.sort()
        merged: list[list[float]] = []
        for span in spans:
            if merged and merged[-1][2] == span[2] and span[0] <= merged[-1][1] + tolerance:
                merged[-1][1] = max(merged[-1][1], span[1])
            else:
                merged.append(span)
        at = [k for k, span in enumerate(merged) if span[2] == own and abs(span[1]) <= tolerance]
        if not at:
            continue
        merged[at[0]][1] = 0.0
        if not ahead:
            merged = merged[: at[0] + 1]
        other = -1 if faced is None else int(faced[index])
        found[index] = _keyed(merged, own, other, body, int(beneath[index]), plan, devices)
    return found


def _keyed(
    spans: list[list[float]],
    own: int,
    faced: int,
    body: Index,
    handle: int,
    plan: FringePlan,
    devices: int,
) -> _Cutting:
    """Cut of ``spans`` with conductors and bodies numbered by first appearance after
    the wall's own, rounded to the grid of ``plan``, equation (160)."""
    step, reach = plan.step, plan.reach
    labels = {own: 0}
    for span in spans:
        labels.setdefault(int(span[2]), len(labels))
    conductors = list(labels)
    owned = [int(body[k]) for k in conductors]
    if handle >= 0:
        owned.append(int(plan.handle_body[handle]))
        conductors.append(devices + handle)
    bodies: dict[int, int] = {}
    for b in owned:
        bodies.setdefault(b, len(bodies))
    low = np.array([span[0] for span in spans])
    high = np.array([span[1] for span in spans])
    rounded_low, rounded_high = np.round(low / step) * step, np.round(high / step) * step
    if np.any(rounded_high <= rounded_low) or np.any(rounded_low[1:] <= rounded_high[:-1]):
        # Silicon or a gap narrower than the grid: this cut keeps its own geometry.
        step = KEY_GRID
        rounded_low, rounded_high = low.copy(), high.copy()
    position = np.zeros(len(bodies))
    offset = np.zeros(len(bodies))
    nearest = np.full(len(bodies), np.inf)
    for k, span in enumerate(spans):
        slot = bodies[int(body[int(span[2])])]
        for actual, rounded in ((low[k], rounded_low[k]), (high[k], rounded_high[k])):
            if abs(abs(actual) - reach) > 1e-6 * reach and abs(actual) < nearest[slot]:
                nearest[slot], position[slot], offset[slot] = abs(actual), actual, actual - rounded
    label = np.array([labels[int(span[2])] for span in spans], dtype=np.intp)
    body_label = np.array([bodies[b] for b in owned], dtype=np.intp)
    cut = Cut(rounded_low, rounded_high, label, body_label, handle >= 0, faced >= 0)
    key = b"".join(
        [
            np.round(rounded_low / step).astype(np.int64).tobytes(),
            np.round(rounded_high / step).astype(np.int64).tobytes(),
            label.tobytes(),
            body_label.tobytes(),
            bytes([handle >= 0, faced >= 0]),
        ]
    )
    faced_body = bodies.get(int(body[faced]), 0) if faced >= 0 and faced in labels else 0
    return _Cutting(
        key,
        cut,
        np.array(conductors, dtype=np.intp),
        np.array(list(bodies), dtype=np.intp),
        position,
        offset,
        faced_body,
    )


def _assemble(
    plan: FringePlan,
    items: _Items,
    pivots: Array,
    centers: Array,
    motion: Array,
    total: Array,
    slope: Array,
    bend: Array,
) -> None:
    """Add each weight times the section of its cut expanded to where its bodies are,
    equations (162) and (163), to the counts on the cores of the items' conductors.

    The variables of an expansion are, for every body of a cut but the wall's, how far
    its nearest wall lies from the rounding along the normal, measured from the wall,
    and for every body its rise, equation (160). Parameters of one item are six per
    block: four blocks for the bodies of its weight and of the exact gap, whose first
    three parameters they use, then one per body of its cut.
    """
    count = len(items.cuttings)
    if not count:
        return
    sections = [plan.sections[cutting.key] for cutting in items.cuttings]
    labels = max(len(section.value) for section in sections)
    held = max(len(cutting.bodies) for cutting in items.cuttings)
    blocks = 4 + held
    params = 6 * blocks
    value: Array = np.zeros((count, labels))
    first: Array = np.zeros((count, labels, held, 2))
    second: Array = np.zeros((count, labels, held, 2, held, 2))
    position: Array = np.zeros((count, held))
    offset: Array = np.zeros((count, held))
    faced: Array = np.zeros(count, dtype=np.intp)
    column = np.full((count, labels), -1, dtype=np.intp)
    member = np.full((count, blocks), -1, dtype=np.intp)
    member[:, :4] = items.slots
    for k, (cutting, section) in enumerate(zip(items.cuttings, sections, strict=True)):
        n, m = len(section.value), len(cutting.bodies)
        value[k, :n] = section.value
        first[k, :n, :m] = section.gradient
        second[k, :n, :m, :, :m] = section.hessian
        position[k, :m] = cutting.position
        offset[k, :m] = cutting.offset
        faced[k] = cutting.faced
        column[k, :n] = cutting.conductors
        member[k, 4 : 4 + m] = cutting.bodies
    present = member[:, 4:] >= 0
    held_body = np.where(present, member[:, 4:], 0)
    # Every body rises at where the cut stands on its wall at rest, as a plate does at
    # its reference point, equations (46) and (160).
    spot = items.start + items.middle.value[:, None] * items.tangent
    lever = spot[:, None, :] - centers[held_body]
    shift: Array = np.zeros((count, held, 2))
    shift[:, 1:, 0] = offset[:, 1:]
    shift[:, :, 1] = (
        motion[held_body, 3]
        + motion[held_body, 4] * lever[:, :, 1]
        - motion[held_body, 5] * lever[:, :, 0]
    )
    shift[~present] = 0.0

    # Rates of the variables in the item's parameters, and their own curvature.
    normal = items.normal
    rate: Array = np.zeros((count, held, 2, blocks, 6))
    curve: Array = np.zeros((count, held, 2, params, params))
    # A wall moves along the cut by where it crosses the line of the cut: turning a
    # body bends that crossing by its arm along the normal, and turning the wall the
    # cut stands on also swings the line, equation (160).
    wall_arm = items.point - pivots[held_body[:, 0]]
    for j in range(1, held):
        their = items.point + normal * position[:, j, None] - pivots[held_body[:, j]]
        for block, arm, sign, bow in (
            (4 + j, their, 1.0, np.einsum("ij,ij->i", normal, their)),
            (4, wall_arm, -1.0, np.einsum("ij,ij->i", normal, wall_arm) + position[:, j]),
        ):
            rate[:, j, 0, block, 0] += sign * normal[:, 0]
            rate[:, j, 0, block, 1] += sign * normal[:, 1]
            rate[:, j, 0, block, 2] += sign * (normal[:, 1] * arm[:, 0] - normal[:, 0] * arm[:, 1])
            where = 6 * block + 2
            curve[:, j, 0, where, where] += bow
    if items.distance is not None:
        gap_rate = np.zeros((count, blocks, 6))
        gap_rate[:, :4, :3] = items.distance.gradient.reshape(count, 4, 3)
        gap_curve = np.zeros((count, blocks, 6, blocks, 6))
        gap_curve[:, :4, :3, :4, :3] = items.distance.hessian.reshape(count, 4, 3, 4, 3)
        on = np.flatnonzero(faced > 0)
        rate[on, faced[on], 0] = gap_rate[on]
        curve[on, faced[on], 0] = gap_curve[on].reshape(len(on), params, params)
    rate[:, :, 1, 4:, 3] = np.eye(held)[None]
    rate[:, :, 1, 4:, 4] = np.eye(held)[None] * lever[:, :, 1, None]
    rate[:, :, 1, 4:, 5] = -np.eye(held)[None] * lever[:, :, 0, None]
    # Where the cut stands moves with the ends of its piece, and turns the lever.
    tip, tilt_x, tilt_y = items.tangent, motion[held_body, 4], motion[held_body, 5]
    along = items.middle.gradient.reshape(count, 4, 3)
    lean = tilt_x * tip[:, None, 1] - tilt_y * tip[:, None, 0]
    rate[:, :, 1, :4, :3] += lean[:, :, None, None] * along[:, None]
    flat_along: Array = np.zeros((count, blocks, 6))
    flat_along[:, :4, :3] = along
    flat_along = flat_along.reshape(count, params)
    for j in range(held):
        for param, arm in ((4, tip[:, 1]), (5, -tip[:, 0])):
            where = 6 * (4 + j) + param
            curve[:, j, 1, where, :] += arm[:, None] * flat_along
            curve[:, j, 1, :, where] += arm[:, None] * flat_along
    bowed: Array = np.zeros((count, blocks, 6, blocks, 6))
    bowed[:, :4, :3, :4, :3] = items.middle.hessian.reshape(count, 4, 3, 4, 3)
    curve[:, :, 1] += lean[:, :, None, None] * bowed.reshape(count, 1, params, params)
    rate[~present] = 0.0
    curve[~present] = 0.0
    variables = 2 * held
    rate = rate.reshape(count, variables, params)
    curve = curve.reshape(count, variables, params * params)
    second = second.reshape(count, labels, variables, variables)
    shift = shift.reshape(count, variables, 1)
    first = first.reshape(count, labels, variables)

    pulled = first + (second @ shift[:, None])[..., 0]
    excess = value + 0.5 * ((first + pulled) @ shift)[..., 0]
    tilt = pulled @ rate
    swing = np.swapaxes(rate, 1, 2)[:, None] @ (second @ rate[:, None])
    swing += (pulled @ curve).reshape(count, labels, params, params)

    run: Array = np.zeros((count, blocks, 6))
    run[:, :4, :3] = items.weight.gradient.reshape(count, 4, 3)
    run_bend: Array = np.zeros((count, blocks, 6, blocks, 6))
    run_bend[:, :4, :3, :4, :3] = items.weight.hessian.reshape(count, 4, 3, 4, 3)
    run = run.reshape(count, params)
    run_bend = run_bend.reshape(count, params, params)
    w = items.weight.value
    gains = w[:, None] * excess
    slopes = run[:, None, :] * excess[:, :, None] + w[:, None, None] * tilt
    cross = run[:, None, :, None] * tilt[:, :, None, :]
    bends = (
        run_bend[:, None] * excess[:, :, None, None]
        + w[:, None, None, None] * swing
        + cross
        + np.swapaxes(cross, 2, 3)
    )
    spare = slope.shape[2] - 1
    target = np.where(member[:, :, None] >= 0, 6 * member[:, :, None] + np.arange(6), spare)
    target = target.reshape(count, params)
    real = column >= 0
    row = np.broadcast_to(items.owner[:, None], column.shape)[real]
    col = column[real]
    places = np.broadcast_to(target[:, None, :], (*column.shape, params))[real]
    np.add.at(total, (row, col), gains[real])
    groups, inverse = np.unique(np.column_stack([row, col, places]), axis=0, return_inverse=True)
    inverse = inverse.ravel()
    order = np.argsort(inverse, kind="stable")
    starts = np.flatnonzero(np.diff(np.concatenate([[-1], inverse[order]])))
    summed_slope = np.add.reduceat(slopes[real][order], starts, axis=0)
    summed_bend = np.add.reduceat(bends[real][order], starts, axis=0)
    for g, group in enumerate(groups):
        r, c, where = int(group[0]), int(group[1]), group[2:]
        np.add.at(slope[r, c], where, summed_slope[g])
        np.add.at(bend[r, c], (where[:, None], where[None, :]), summed_bend[g])


def _symmetric(
    total: Array, slope: Array, bend: Array, devices: int, bodies: int
) -> FringeCapacitance:
    """Counts on the cores of device electrodes as a symmetric matrix with rows that sum
    to zero: entries between two device electrodes are the mean of both counts, those
    to the handle the count on the device electrode alone."""
    count = len(total)
    weight = np.ones((count, count))
    weight[:devices, :devices] = 2.0
    diagonal = np.arange(count)
    parts = []
    for part in (total, slope, bend):
        both = (part + np.swapaxes(part, 0, 1)) / weight.reshape(
            count, count, *([1] * (part.ndim - 2))
        )
        both[diagonal, diagonal] = 0.0
        both[diagonal, diagonal] = -both.sum(axis=1)
        parts.append(both)
    return FringeCapacitance(
        parts[0],
        parts[1].reshape(count, count, bodies, 6),
        parts[2].reshape(count, count, bodies, 6, bodies, 6),
    )

"""Capacitance between conductor outlines in the layout plane."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import scipy.linalg
from numpy.typing import ArrayLike, NDArray
from shapely import MultiPolygon, Polygon, STRtree

from connflex_core.electrostatics import VACUUM_PERMITTIVITY
from connflex_core.electrostatics.kernel import pair_derivatives, pair_integrals
from connflex_core.electrostatics.outlines import Mesh, mesh_outlines
from connflex_core.geometry.heal import as_multipolygon
from connflex_core.params import ElectrostaticParams
from connflex_core.progress import Progress, stage

Array = NDArray[np.float64]

PAIRS_PER_PASS = 250_000
DERIVATIVE_PAIRS_PER_PASS = 40_000


@dataclass(frozen=True, eq=False)
class SidewallCapacitance:
    """Capacitance between conductor outlines in the layout plane, per meter of thickness.

    Attributes:
        capacitance: n by n matrix in farads per meter: charge per meter of thickness on
            one conductor per volt on another. Symmetric, and every row sums to zero.
        gradient: n by n by m by 3 derivatives of ``capacitance`` for travel of each of m
            rigid bodies along x and y, in farads per square meter, and for its rotation
            about z through its center, in farads per meter and radian.
        hessian: n by n by m by 3 by m by 3 second derivatives of ``capacitance`` for one
            such motion of each of two bodies. The point a body rotates about travels
            with it.
        centers: m by 2 points about which the bodies rotate, in meters.
        charge: N by n charges on the elements of ``mesh``, in coulombs per meter of
            thickness, per volt on each conductor with the others at zero.
        mesh: Boundary elements on the outlines.
    """

    capacitance: Array
    gradient: Array
    hessian: Array
    centers: Array
    charge: Array
    mesh: Mesh


def sidewall_capacitance(
    conductors: Sequence[Polygon | MultiPolygon],
    params: ElectrostaticParams,
    *,
    bodies: Sequence[int] | None = None,
    centers: ArrayLike | None = None,
    progress: Progress | None = None,
) -> SidewallCapacitance:
    """Capacitance between sidewalls of conductors drawn in the layout plane.

    Each entry of ``conductors`` is one equipotential region, in meters; regions must
    not overlap or touch. The field is two-dimensional, the set of conductors carries no
    net charge, and potential far away floats, equations (61) to (71). ``params`` must
    carry ``element_size``.

    ``bodies`` gives the rigid body of each conductor as an index from zero; conductors
    of one body move together. By default every conductor is its own body. ``centers``
    gives the point each body rotates about, by default the centroid of its conductors.
    Gradient and Hessian are the exact derivatives of the computed capacitance for rigid
    motion of the bodies, equations (75) to (90). Stages mesh, assemble, solve,
    derivatives, and curvature log through the ``connflex_core`` logger and report to
    ``progress``.
    """
    if params.element_size is None:
        raise ValueError("element_size is required for sidewall capacitance")
    regions = _disjoint(conductors)
    body = np.arange(len(regions)) if bodies is None else np.asarray(bodies, dtype=np.intp)
    rigid = int(np.max(body, initial=-1)) + 1
    if body.shape != (len(regions),) or not np.array_equal(np.unique(body), np.arange(rigid)):
        raise ValueError("bodies must give every conductor an index, using each index from zero up")
    if centers is None:
        members = [
            [
                polygon
                for region, b in zip(regions, body, strict=True)
                if b == k
                for polygon in region.geoms
            ]
            for k in range(rigid)
        ]
        pivots = np.array([MultiPolygon(group).centroid.coords[0] for group in members])
    else:
        pivots = np.asarray(centers, dtype=float)
        if pivots.shape != (rigid, 2):
            raise ValueError(f"centers must hold one (x, y) point for each of {rigid} bodies")
    with stage("mesh", progress) as current:
        mesh = mesh_outlines(regions, params.element_size, params.corner_grading)
        count = len(mesh.conductor)
        current.count(elements=count)
    if count > params.max_elements:
        raise ValueError(
            f"outlines need {count} boundary elements, more than max_elements = "
            f"{params.max_elements}: raise element_size or max_elements, or pass fewer conductors"
        )

    corners = np.concatenate([mesh.start, mesh.end])
    low, high = corners.min(axis=0), corners.max(axis=0)
    scale = 2 * float(np.hypot(*(high - low)))
    start, end = (
        (points - (low + high) / 2) / scale @ (1, 1j) for points in (mesh.start, mesh.end)
    )
    length = np.abs(end - start)

    potential = np.empty((count, count))
    rows = max(1, PAIRS_PER_PASS // count)
    with stage("assemble", progress, total=count) as current:
        for first in range(0, count, rows):
            block = slice(first, first + rows)
            integral = pair_integrals(start[block, None], end[block, None], start, end)
            potential[block] = -integral / (2 * np.pi * length[block, None] * length)
            current.advance(len(length[block]))

    with stage("solve", progress):
        try:
            factor = scipy.linalg.cho_factor(potential, lower=True, overwrite_a=True)
        except scipy.linalg.LinAlgError as exc:
            raise ValueError("outlines cross or elements coincide") from exc
        on = (mesh.conductor[:, None] == np.arange(len(regions))).astype(float)
        solved = scipy.linalg.cho_solve(factor, np.column_stack([on, np.ones(count)]))
        driven, uniform = solved[:, :-1], solved[:, -1]
        charge = driven - np.outer(uniform, uniform @ on) / uniform.sum()
        capacitance = on.T @ charge

    density = charge / length[:, None]
    carried = [np.flatnonzero(body[mesh.conductor] == k) for k in range(rigid)]
    pivot = (pivots - (low + high) / 2) / scale @ (1, 1j)
    shape = (len(regions), len(regions), 3 * rigid)
    gradient = np.zeros(shape)
    hessian = np.zeros((*shape, 3 * rigid))
    swing = np.zeros((count, len(regions), 3 * rigid))
    pairs = [(k, other) for k in range(rigid) for other in range(k + 1, rigid)]
    with stage("derivatives", progress, total=len(pairs)) as current:
        for k, other in pairs:
            mine, theirs = carried[k], carried[other]
            freedom = np.r_[3 * k : 3 * k + 3, 3 * other : 3 * other + 3]
            rows = max(1, DERIVATIVE_PAIRS_PER_PASS // theirs.size)
            for first in range(0, mine.size, rows):
                part = mine[first : first + rows]
                slopes, bends = pair_derivatives(
                    start[part, None],
                    end[part, None],
                    start[theirs],
                    end[theirs],
                    pivot[k],
                    pivot[other],
                )
                for p, slope in enumerate(slopes):
                    onto = slope @ density[theirs]
                    gradient[:, :, freedom[p]] += density[part].T @ onto
                    swing[part, :, freedom[p]] += onto / length[part, None]
                    swing[theirs, :, freedom[p]] += slope.T @ density[part] / length[theirs, None]
                    for q in range(6):
                        hessian[:, :, freedom[p], freedom[q]] += (
                            density[part].T @ bends[p, q] @ density[theirs]
                        )
            current.advance()
    gradient = (gradient + gradient.transpose(1, 0, 2)) / (2 * np.pi)
    hessian = (hessian + hessian.transpose(1, 0, 2, 3)) / (2 * np.pi)

    with stage("curvature", progress):
        swing /= -2 * np.pi
        flat = swing.reshape(count, -1)
        settled = scipy.linalg.cho_solve(factor, flat)
        settled -= np.outer(uniform, uniform @ flat) / uniform.sum()
        through = (flat.T @ settled).reshape(*shape[1:], *shape[1:]).transpose(0, 2, 1, 3)
        hessian += through + through.transpose(0, 1, 3, 2)

    unit = np.tile([1 / scale, 1 / scale, 1.0], rigid)
    gradient *= unit
    hessian *= unit[:, None] * unit
    permittivity = VACUUM_PERMITTIVITY * params.ambient_permittivity
    return SidewallCapacitance(
        capacitance=permittivity * (capacitance + capacitance.T) / 2,
        gradient=permittivity * gradient.reshape(len(regions), len(regions), rigid, 3),
        hessian=permittivity * hessian.reshape(len(regions), len(regions), rigid, 3, rigid, 3),
        centers=pivots,
        charge=permittivity * charge,
        mesh=mesh,
    )


def _disjoint(conductors: Sequence[Polygon | MultiPolygon]) -> list[MultiPolygon]:
    regions = []
    for index, conductor in enumerate(conductors):
        region = as_multipolygon(conductor)
        if region.is_empty or not conductor.is_valid:
            raise ValueError(f"conductor {index} is empty or not a valid polygon")
        regions.append(region)
    polygons = [polygon for region in regions for polygon in region.geoms]
    owner = np.array([index for index, region in enumerate(regions) for _ in region.geoms])
    first, second = STRtree(polygons).query(polygons, predicate="intersects")
    clash = owner[first] < owner[second]
    if clash.any():
        a, b = owner[first][clash][0], owner[second][clash][0]
        raise ValueError(f"conductors {a} and {b} overlap or touch")
    return regions

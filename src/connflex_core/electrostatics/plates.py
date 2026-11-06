"""Plates between device and handle electrodes."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import shapely
from numpy.typing import ArrayLike, NDArray
from shapely import MultiPolygon, Polygon, STRtree
from shapely.geometry.base import BaseGeometry

from connflex_core.electrostatics import VACUUM_PERMITTIVITY
from connflex_core.geometry.heal import as_multipolygon
from connflex_core.geometry.release import openings_within
from connflex_core.mechanics.inertia import AreaMoments, area_moments
from connflex_core.params import ElectrostaticParams
from connflex_core.progress import Progress, stage

if TYPE_CHECKING:
    from connflex_core.model import Structure

Array = NDArray[np.float64]


@dataclass(frozen=True, eq=False)
class Plate:
    """Capacitor between a device electrode and the handle electrode beneath it.

    Attributes:
        index: Position in the tuple returned by ``Structure.plates``.
        electrodes: Indices of the device electrode and the handle electrode.
        bodies: Indices of the bodies those electrodes move with.
        oxide: Capacitance through surviving buried oxide, in farads. Motion does not
            change it.
        gap: Height of the released gap at rest, in meters: thickness of the box layer.
        permittivity: Permittivity of the medium in the released gap, in farads per meter.
        origin: Reference point (x, y, z) in meters: center of the bounding box of the
            overlap, at mid-height of the gap.
        extent: Half the width of that bounding box along x and along y.
        moments: 3 by 3 integrals over the released overlap of products of 1, x, and y
            measured from ``origin``, equation (48).
        edges: 3 by 3 integrals along the handle outline under device silicon,
            equation (52): one row for travel along x, along y, and about z, one column
            for weight 1, x, and y measured from ``origin``.
        centers: Centers of mass of the two bodies, one row each.
    """

    index: int
    electrodes: tuple[int, int]
    bodies: tuple[int, int]
    oxide: float
    gap: float
    permittivity: float
    origin: Array
    extent: Array
    moments: Array
    edges: Array
    centers: Array

    @property
    def released(self) -> float:
        """Capacitance across the released gap at rest, in farads."""
        return float(self.permittivity * self.moments[0, 0] / self.gap)

    def lowest(self, device: ArrayLike | None = None, handle: ArrayLike | None = None) -> float:
        """Smallest gap over the bounding box of the overlap, in meters, equation (126).

        ``device`` and ``handle`` are as for ``evaluate``. Zero or less means that a
        corner of the box has reached the other plate.
        """
        local = self._local(device, handle)
        return float(self.gap + local[3] - np.abs(local[4:]) @ self.extent)

    def _local(self, device: ArrayLike | None, handle: ArrayLike | None) -> Array:
        """Travel, rise, and slopes at ``origin`` for displacements of the two bodies."""
        displaced = [np.zeros(6) if q is None else np.asarray(q, float) for q in (device, handle)]
        local: Array = self._jacobian() @ np.concatenate(displaced)
        return local

    def _jacobian(self) -> Array:
        return np.hstack(
            [_local(self.origin, self.centers[0]), -_local(self.origin, self.centers[1])]
        )

    def evaluate(
        self, device: ArrayLike | None = None, handle: ArrayLike | None = None
    ) -> tuple[float, Array, Array]:
        """Capacitance in farads with its gradient and Hessian, equations (55) to (58).

        ``device`` and ``handle`` are the six degrees of freedom of equation (12) of the
        body of each electrode, zero when omitted. The gradient has twelve entries and the
        Hessian is 12 by 12: degrees of freedom of the device electrode's body come first.
        """
        jacobian = self._jacobian()
        local = self._local(device, handle)
        travel, rise, slope = local[:3], local[3], local[4:]
        gap = self.gap + rise
        if gap <= 0:
            raise ValueError(f"gap between plates closes: {gap} m at the reference point")
        area = self.moments[0, 0] - travel @ self.edges[:, 0]
        first = self.moments[0, 1:] - travel @ self.edges[:, 1:]
        second = self.moments[1:, 1:]
        linear, quadratic = slope @ first, slope @ second @ slope
        value = area / gap - linear / gap**2 + quadratic / gap**3
        gradient = np.empty(6)
        gradient[:3] = -self.edges[:, 0] / gap + self.edges[:, 1:] @ slope / gap**2
        gradient[3] = -area / gap**2 + 2 * linear / gap**3 - 3 * quadratic / gap**4
        gradient[4:] = -first / gap**2 + 2 * second @ slope / gap**3
        hessian = np.zeros((6, 6))
        hessian[:3, 3] = self.edges[:, 0] / gap**2 - 2 * self.edges[:, 1:] @ slope / gap**3
        hessian[:3, 4:] = self.edges[:, 1:] / gap**2
        hessian[3, 3] = 2 * area / gap**3 - 6 * linear / gap**4 + 12 * quadratic / gap**5
        hessian[3, 4:] = 2 * first / gap**3 - 6 * second @ slope / gap**4
        hessian[4:, 4:] = 2 * second / gap**3
        hessian = np.triu(hessian) + np.triu(hessian, 1).T
        return (
            float(self.oxide + self.permittivity * value),
            self.permittivity * jacobian.T @ gradient,
            self.permittivity * jacobian.T @ hessian @ jacobian,
        )


def find_plates(
    structure: Structure,
    params: ElectrostaticParams | None = None,
    *,
    progress: Progress | None = None,
) -> tuple[Plate, ...]:
    """Capacitors between device electrodes and handle electrodes beneath them.

    One plate per pair of a device and a handle electrode that overlap, ordered by device
    electrode and then by handle electrode. Its overlap is split by equation (42) into
    oxide and released gap; release openings carry no charge. Pieces of handle outline
    under device silicon give the change of overlap with in-plane travel. Stage plates
    logs through the ``connflex_core`` logger and reports to ``progress`` once per device
    electrode.
    """
    params = params or ElectrostaticParams()
    stack = structure.stack
    relative = stack.box.material.relative_permittivity
    if relative is None:
        raise ValueError("box material has no relative_permittivity")
    gap = stack.box.thickness
    height = (stack.z_extents[0][1] + stack.z_extents[2][0]) / 2
    device = [e for e in structure.electrodes if e.layer == "device"]
    handle = [e for e in structure.electrodes if e.layer == "handle"]
    plates: list[Plate] = []
    with stage("plates", progress, total=len(device)) as current:
        below = STRtree([electrode.region for electrode in handle])
        oxide = [polygon for body in structure.bodies for polygon in body.box.geoms]
        openings = list(structure.openings.geoms)
        oxide_tree, opening_tree = STRtree(oxide), STRtree(openings)
        for above in device:
            silicon = shapely.union_all(shapely.get_parts(above.region))
            for j in sorted(below.query(silicon, predicate="intersects").tolist()):
                under = handle[j]
                overlap = as_multipolygon(silicon.intersection(under.region))
                if overlap.is_empty:
                    continue
                bonded = area_moments(MultiPolygon(openings_within(oxide_tree, oxide, overlap)))
                opened = area_moments(
                    MultiPolygon(openings_within(opening_tree, openings, overlap))
                )
                west, south, east, north = overlap.bounds
                origin = np.array([(west + east) / 2, (south + north) / 2, height])
                outline = _crossings(under.region, silicon, openings, opening_tree)
                plates.append(
                    Plate(
                        index=len(plates),
                        electrodes=(above.index, under.index),
                        bodies=(above.body, under.body),
                        oxide=VACUUM_PERMITTIVITY * relative * bonded.area / gap,
                        gap=gap,
                        permittivity=VACUUM_PERMITTIVITY * params.ambient_permittivity,
                        origin=origin,
                        extent=np.array([(east - west) / 2, (north - south) / 2]),
                        moments=_about(area_moments(overlap) - bonded - opened, origin),
                        edges=_edge_integrals(*outline, origin),
                        centers=np.array(
                            [structure.bodies[i].inertial.center for i in (above.body, under.body)]
                        ),
                    )
                )
            current.advance()
        current.count(plates=len(plates))
    return tuple(plates)


def _about(moments: AreaMoments, origin: Array) -> Array:
    """Moment matrix of equation (48) about ``origin``, by the shift of equation (11c)."""
    x, y = (float(value) for value in origin[:2])
    m = moments
    sx, sy = m.x - x * m.area, m.y - y * m.area
    sxx = m.xx - 2 * x * m.x + x * x * m.area
    syy = m.yy - 2 * y * m.y + y * y * m.area
    sxy = m.xy - x * m.y - y * m.x + x * y * m.area
    return np.array([[m.area, sx, sy], [sx, sxx, sxy], [sy, sxy, syy]])


def _crossings(
    under: MultiPolygon,
    silicon: BaseGeometry,
    openings: Sequence[Polygon],
    opening_tree: STRtree,
) -> tuple[Array, Array, Array]:
    """Straight pieces of the outline of ``under`` beneath ``silicon`` less ``openings``, (50).

    Returns start points, end points, and the unit normal out of ``under`` of each piece.
    """
    rings, polygon_of = shapely.get_rings(shapely.get_parts(under), return_index=True)
    coords, ring_of = shapely.get_coordinates(rings, return_index=True)
    edge = ring_of[:-1] == ring_of[1:]
    start, end, ring = coords[:-1][edge], coords[1:][edge], ring_of[:-1][edge]
    cross = start[:, 0] * end[:, 1] - end[:, 0] * start[:, 1]
    exterior = np.ones(rings.size, dtype=bool)
    exterior[1:] = polygon_of[1:] != polygon_of[:-1]
    counterclockwise = np.bincount(ring, weights=cross, minlength=rings.size) > 0
    outward = np.where(exterior == counterclockwise, 1.0, -1.0)[ring]
    along = end - start
    length = np.hypot(along[:, 0], along[:, 1])
    kept = length > 0
    normals = outward[kept, None] * np.stack([along[kept, 1], -along[kept, 0]], axis=1)
    normals /= length[kept, None]

    lines = np.asarray(shapely.linestrings(np.stack([start[kept], end[kept]], axis=1)))
    shapely.prepare(silicon)
    inside = np.flatnonzero(shapely.intersects(lines, silicon))
    pieces = np.asarray(shapely.intersection(lines[inside], silicon))
    pierced, hit = opening_tree.query(pieces, predicate="intersects")
    for i in np.unique(pierced):
        holes = shapely.union_all([openings[j] for j in hit[pierced == i]])
        pieces[i] = shapely.difference(pieces[i], holes)
    parts, piece_of = shapely.get_parts(pieces, return_index=True)
    parts, part_of = shapely.get_parts(parts, return_index=True)
    points, line_of = shapely.get_coordinates(parts, return_index=True)
    joined = line_of[:-1] == line_of[1:]
    owner = inside[piece_of[part_of[line_of[:-1][joined]]]]
    return points[:-1][joined], points[1:][joined], normals[owner]


def _edge_integrals(start: Array, end: Array, normal: Array, origin: Array) -> Array:
    """Integrals of equation (52) over straight pieces, in the closed form of (54)."""
    from_origin = [start - origin[:2], end - origin[:2]]
    length = np.hypot(*(end - start).T)
    flow = [
        np.stack([normal[:, 0], normal[:, 1], normal[:, 1] * r[:, 0] - normal[:, 0] * r[:, 1]])
        for r in from_origin
    ]
    weight = [np.stack([np.ones(len(r)), r[:, 0], r[:, 1]]) for r in from_origin]
    edges: Array = (
        np.einsum("m,km,jm->kj", length, 2 * flow[0] + flow[1], weight[0])
        + np.einsum("m,km,jm->kj", length, flow[0] + 2 * flow[1], weight[1])
    ) / 6
    return edges


def _local(origin: Array, center: Array) -> Array:
    """Jacobian of travel, rise, and slopes at ``origin`` for a body with ``center``, (46)."""
    dx, dy, dz = origin - center
    return np.array(
        [
            [1.0, 0.0, 0.0, 0.0, dz, -dy],
            [0.0, 1.0, 0.0, -dz, 0.0, dx],
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
            [0.0, 0.0, 1.0, dy, -dx, 0.0],
            [0.0, 0.0, 0.0, 0.0, -1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
        ]
    )

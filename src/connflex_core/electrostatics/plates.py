"""Plates between device and handle electrodes."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import shapely
from numpy.typing import ArrayLike, NDArray
from shapely import MultiPolygon, STRtree

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
        origin: Reference point (x, y) in meters: center of the bounding box of the overlap.
        moments: 3 by 3 integrals over the released overlap of products of 1, x, and y
            measured from ``origin``, equation (48).
        centers: Centers of mass of the two bodies, one row each.
    """

    index: int
    electrodes: tuple[int, int]
    bodies: tuple[int, int]
    oxide: float
    gap: float
    permittivity: float
    origin: Array
    moments: Array
    centers: Array

    @property
    def released(self) -> float:
        """Capacitance across the released gap at rest, in farads."""
        return float(self.permittivity * self.moments[0, 0] / self.gap)

    def evaluate(
        self, device: ArrayLike | None = None, handle: ArrayLike | None = None
    ) -> tuple[float, Array, Array]:
        """Capacitance in farads with its gradient and Hessian, equations (49) to (52).

        ``device`` and ``handle`` are the six degrees of freedom of equation (12) of the
        body of each electrode, zero when omitted. The gradient has twelve entries and the
        Hessian is 12 by 12: degrees of freedom of the device electrode's body come first.
        """
        jacobian = np.hstack(
            [_local(self.origin, self.centers[0]), -_local(self.origin, self.centers[1])]
        )
        displaced = [np.zeros(6) if q is None else np.asarray(q, float) for q in (device, handle)]
        rise, *slopes = jacobian @ np.concatenate(displaced)
        gap = self.gap + rise
        if gap <= 0:
            raise ValueError(f"gap between plates closes: {gap} m at the reference point")
        slope = np.array(slopes)
        area, first, second = self.moments[0, 0], self.moments[0, 1:], self.moments[1:, 1:]
        linear, quadratic = slope @ first, slope @ second @ slope
        value = area / gap - linear / gap**2 + quadratic / gap**3
        gradient = np.empty(3)
        gradient[0] = -area / gap**2 + 2 * linear / gap**3 - 3 * quadratic / gap**4
        gradient[1:] = -first / gap**2 + 2 * second @ slope / gap**3
        hessian = np.empty((3, 3))
        hessian[0, 0] = 2 * area / gap**3 - 6 * linear / gap**4 + 12 * quadratic / gap**5
        hessian[0, 1:] = hessian[1:, 0] = 2 * first / gap**3 - 6 * second @ slope / gap**4
        hessian[1:, 1:] = 2 * second / gap**3
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
    oxide and released gap; release openings carry no charge. Stage plates logs through
    the ``connflex_core`` logger and reports to ``progress`` once per device electrode.
    """
    params = params or ElectrostaticParams()
    stack = structure.stack
    if stack is None:
        raise ValueError("structure has no stack")
    relative = stack.box.material.relative_permittivity
    if relative is None:
        raise ValueError("box material has no relative_permittivity")
    gap = stack.box.thickness
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
                origin = np.array([(west + east) / 2, (south + north) / 2])
                plates.append(
                    Plate(
                        index=len(plates),
                        electrodes=(above.index, under.index),
                        bodies=(above.body, under.body),
                        oxide=VACUUM_PERMITTIVITY * relative * bonded.area / gap,
                        gap=gap,
                        permittivity=VACUUM_PERMITTIVITY * params.ambient_permittivity,
                        origin=origin,
                        moments=_about(area_moments(overlap) - bonded - opened, origin),
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
    x, y = (float(value) for value in origin)
    m = moments
    sx, sy = m.x - x * m.area, m.y - y * m.area
    sxx = m.xx - 2 * x * m.x + x * x * m.area
    syy = m.yy - 2 * y * m.y + y * y * m.area
    sxy = m.xy - x * m.y - y * m.x + x * y * m.area
    return np.array([[m.area, sx, sy], [sx, sxx, sxy], [sy, sxy, syy]])


def _local(origin: Array, center: Array) -> Array:
    """Jacobian of gap rise and slopes at ``origin`` in degrees of freedom of a body, (46)."""
    dx, dy = origin - center[:2]
    return np.array(
        [
            [0.0, 0.0, 1.0, dy, -dx, 0.0],
            [0.0, 0.0, 0.0, 0.0, -1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
        ]
    )

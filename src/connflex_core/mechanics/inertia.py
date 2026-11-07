"""Mass properties of bodies."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import shapely
from numpy.typing import NDArray
from shapely import MultiPolygon, Polygon

from connflex_core.stack import Stack


@dataclass(frozen=True)
class AreaMoments:
    """Integrals of 1, x, y, x^2, y^2, and xy over a region, about the origin.

    Attributes:
        area: Integral of 1.
        x: Integral of x.
        y: Integral of y.
        xx: Integral of x squared.
        yy: Integral of y squared.
        xy: Integral of x times y.
    """

    area: float
    x: float
    y: float
    xx: float
    yy: float
    xy: float

    def about(self, x: float, y: float) -> AreaMoments:
        """These moments about the point (x, y), by the shift of equation (11c)."""
        return AreaMoments(
            self.area,
            self.x - x * self.area,
            self.y - y * self.area,
            self.xx - 2 * x * self.x + x * x * self.area,
            self.yy - 2 * y * self.y + y * y * self.area,
            self.xy - x * self.y - y * self.x + x * y * self.area,
        )

    def __sub__(self, other: AreaMoments) -> AreaMoments:
        return AreaMoments(
            self.area - other.area,
            self.x - other.x,
            self.y - other.y,
            self.xx - other.xx,
            self.yy - other.yy,
            self.xy - other.xy,
        )


@dataclass(frozen=True, eq=False)
class MassProperties:
    """Mass, center of mass, and inertia tensor of a body.

    Attributes:
        mass: Mass in kilograms.
        center: Center of mass (x, y, z) in meters, global frame.
        inertia: 3 by 3 inertia tensor about the center of mass, axes parallel to global
            axes, in kilogram square meters.
    """

    mass: float
    center: NDArray[np.float64]
    inertia: NDArray[np.float64]

    @property
    def matrix(self) -> NDArray[np.float64]:
        """6 by 6 mass matrix in body degrees of freedom, equation (14)."""
        matrix = np.zeros((6, 6))
        matrix[:3, :3] = self.mass * np.eye(3)
        matrix[3:, 3:] = self.inertia
        return matrix


def area_moments(geometry: MultiPolygon) -> AreaMoments:
    """Area moments in closed form over the polygon edges, equations (11a) and (11b).

    Exterior rings add and holes subtract, whichever way they are drawn. The edge sums
    run about the center of the bounding box, where they cancel least, and are moved to
    the origin by the shift of (11c).
    """
    if geometry.is_empty:
        return AreaMoments(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    polygons = shapely.get_parts(geometry)
    rings, owner = shapely.get_rings(polygons, return_index=True)
    coords, ring_of = shapely.get_coordinates(rings, return_index=True)
    first = np.ones(rings.size, dtype=bool)
    first[1:] = owner[1:] != owner[:-1]
    role = np.where(first, 1.0, -1.0)
    shift = (coords.min(axis=0) + coords.max(axis=0)) / 2
    coords = coords - shift

    edge = ring_of[:-1] == ring_of[1:]
    x0, y0 = coords[:-1, 0][edge], coords[:-1, 1][edge]
    x1, y1 = coords[1:, 0][edge], coords[1:, 1][edge]
    ring = ring_of[:-1][edge]
    cross = x0 * y1 - x1 * y0
    terms = np.stack(
        [
            cross / 2,
            (x0 + x1) * cross / 6,
            (y0 + y1) * cross / 6,
            (x0 * x0 + x0 * x1 + x1 * x1) * cross / 12,
            (y0 * y0 + y0 * y1 + y1 * y1) * cross / 12,
            (2 * x0 * y0 + x0 * y1 + x1 * y0 + 2 * x1 * y1) * cross / 24,
        ]
    )
    per_ring = np.stack([np.bincount(ring, weights=row, minlength=rings.size) for row in terms])
    signs = role * np.sign(per_ring[0])
    area, x, y, xx, yy, xy = (float(value) for value in per_ring @ signs)
    sx, sy = (float(value) for value in shift)
    return AreaMoments(area, x, y, xx, yy, xy).about(-sx, -sy)


def mass_properties(
    layers: Sequence[MultiPolygon], stack: Stack, openings: Sequence[Polygon] = ()
) -> MassProperties:
    """Mass properties of a body from its polygons in each stack layer, listed top down.

    Each layer is a uniform extrusion between its top and bottom heights, equations (5)
    to (10). ``openings`` are disjoint release openings lying inside the first layer; their
    moments are subtracted from it, which equals integrating over the perforated polygon
    and costs one pass over the vertices of the openings.
    """
    moments = [area_moments(geometry) for geometry in layers]
    if openings:
        moments[0] = moments[0] - area_moments(MultiPolygon(openings))
    densities = [layer.material.density for layer in stack.layers]
    extents = stack.z_extents
    masses = [
        rho * (top - bottom) * m.area
        for rho, (top, bottom), m in zip(densities, extents, moments, strict=True)
    ]
    mass = sum(masses)
    if mass <= 0:
        raise ValueError("body has no mass")
    cx = (
        sum(
            rho * (top - bottom) * m.x
            for rho, (top, bottom), m in zip(densities, extents, moments, strict=True)
        )
        / mass
    )
    cy = (
        sum(
            rho * (top - bottom) * m.y
            for rho, (top, bottom), m in zip(densities, extents, moments, strict=True)
        )
        / mass
    )
    cz = (
        sum(m_l * (top + bottom) / 2 for m_l, (top, bottom) in zip(masses, extents, strict=True))
        / mass
    )

    inertia = np.zeros((3, 3))
    for rho, (top, bottom), m in zip(densities, extents, moments, strict=True):
        if m.area == 0:
            continue
        shifted = m.about(cx, cy)
        sx, sy, sxx, syy, sxy = shifted.x, shifted.y, shifted.xx, shifted.yy, shifted.xy
        zt, zb = top - cz, bottom - cz
        dz, z2, z3 = zt - zb, (zt * zt - zb * zb) / 2, (zt**3 - zb**3) / 3
        inertia[0, 0] += rho * (dz * syy + z3 * m.area)
        inertia[1, 1] += rho * (dz * sxx + z3 * m.area)
        inertia[2, 2] += rho * dz * (sxx + syy)
        inertia[0, 1] -= rho * dz * sxy
        inertia[0, 2] -= rho * z2 * sx
        inertia[1, 2] -= rho * z2 * sy
    inertia[1, 0], inertia[2, 0], inertia[2, 1] = inertia[0, 1], inertia[0, 2], inertia[1, 2]
    return MassProperties(mass=mass, center=np.array([cx, cy, cz]), inertia=inertia)

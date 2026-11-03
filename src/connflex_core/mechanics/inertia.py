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
    """Area moments by constrained Delaunay triangulation and edge-midpoint rule.

    The midpoint rule integrates polynomials of degree two exactly on each triangle,
    equation (11).
    """
    if geometry.is_empty:
        return AreaMoments(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    triangles = shapely.get_parts(shapely.constrained_delaunay_triangles(geometry))
    vertices = shapely.get_coordinates(triangles).reshape(-1, 4, 2)[:, :3, :]
    v0, v1, v2 = vertices[:, 0], vertices[:, 1], vertices[:, 2]
    cross = (v1[:, 0] - v0[:, 0]) * (v2[:, 1] - v0[:, 1]) - (v2[:, 0] - v0[:, 0]) * (
        v1[:, 1] - v0[:, 1]
    )
    weights = np.tile(np.abs(cross) / 6.0, 3)
    midpoints = np.concatenate([(v0 + v1) / 2, (v1 + v2) / 2, (v2 + v0) / 2])
    x, y = midpoints[:, 0], midpoints[:, 1]
    return AreaMoments(
        area=float(weights.sum()),
        x=float(weights @ x),
        y=float(weights @ y),
        xx=float(weights @ (x * x)),
        yy=float(weights @ (y * y)),
        xy=float(weights @ (x * y)),
    )


def mass_properties(
    layers: Sequence[MultiPolygon], stack: Stack, openings: Sequence[Polygon] = ()
) -> MassProperties:
    """Mass properties of a body from its polygons in each stack layer, listed top down.

    Each layer is a uniform extrusion between its top and bottom heights, equations (5)
    to (10). ``openings`` are disjoint release openings lying inside the first layer; their
    moments are subtracted from it, which equals integrating over the perforated polygon
    without triangulating a polygon with thousands of holes.
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
        sx, sy = m.x - cx * m.area, m.y - cy * m.area
        sxx = m.xx - 2 * cx * m.x + cx * cx * m.area
        syy = m.yy - 2 * cy * m.y + cy * cy * m.area
        sxy = m.xy - cx * m.y - cy * m.x + cx * cy * m.area
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

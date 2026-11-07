"""Rigid motion of points on a body."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

Array = NDArray[np.float64]


def skew(vector: ArrayLike) -> Array:
    """Matrix of the cross product with ``vector`` from the left."""
    x, y, z = np.asarray(vector, dtype=float)
    return np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])


def point_jacobian(offset: ArrayLike) -> Array:
    """Displacement of a point ``offset`` from a body's center per unit translation and
    small rotation of the body, 3 by 6: the identity beside minus the skew of the
    offset, equations (27) and (46)."""
    return np.hstack([np.eye(3), -skew(offset)])

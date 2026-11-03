"""System matrices on body degrees of freedom."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from connflex_core.mechanics.stiffness import effective_stiffness
from connflex_core.model import Structure

DOF_NAMES = ("ux", "uy", "uz", "tx", "ty", "tz")
Dof = tuple[int, str]


@dataclass(frozen=True, eq=False)
class SystemMatrices:
    """Assembled matrices over the movable bodies, equations (34) and (35).

    Attributes:
        mass: 6N by 6N block-diagonal mass matrix.
        stiffness: 6N by 6N stiffness matrix.
        dofs: (body index, degree of freedom) for each row, six per movable body in body
            order: ux, uy, uz, tx, ty, tz.
    """

    mass: NDArray[np.float64]
    stiffness: NDArray[np.float64]
    dofs: tuple[Dof, ...]


def assemble(structure: Structure) -> SystemMatrices:
    """Assemble the mass and stiffness matrices; the anchor's rows and columns are dropped."""
    movable = [body for body in structure.bodies if not body.anchor]
    offset = {body.index: 6 * position for position, body in enumerate(movable)}
    size = 6 * len(movable)
    mass = np.zeros((size, size))
    stiffness = np.zeros((size, size))
    for body in movable:
        row = offset[body.index]
        mass[row : row + 6, row : row + 6] = body.inertial.matrix
    centers = {body.index: body.inertial.center for body in structure.bodies}
    for beam in structure.beams:
        a, b = beam.bodies
        k = effective_stiffness(beam, structure.material, centers[a], centers[b])
        for i, body_i in enumerate((a, b)):
            for j, body_j in enumerate((a, b)):
                if body_i in offset and body_j in offset:
                    ri, rj = offset[body_i], offset[body_j]
                    stiffness[ri : ri + 6, rj : rj + 6] += k[6 * i : 6 * i + 6, 6 * j : 6 * j + 6]
    dofs = tuple((body.index, name) for body in movable for name in DOF_NAMES)
    return SystemMatrices(mass=mass, stiffness=stiffness, dofs=dofs)

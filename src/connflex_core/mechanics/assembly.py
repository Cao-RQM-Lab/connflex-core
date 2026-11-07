"""System matrices on body degrees of freedom."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import networkx as nx
import numpy as np
from numpy.typing import NDArray

from connflex_core.mechanics.stiffness import effective_stiffness
from connflex_core.model import Structure
from connflex_core.progress import Progress, stage

DOF_NAMES = ("ux", "uy", "uz", "tx", "ty", "tz")
Dof = tuple[int, str]


@dataclass(frozen=True, eq=False)
class SystemMatrices:
    """Assembled matrices over movable bodies, equations (34) and (35).

    Attributes:
        mass: 6N by 6N block-diagonal mass matrix.
        stiffness: 6N by 6N stiffness matrix.
        dofs: (body index, degree of freedom) for each row, six per movable body in body
            order: ux, uy, uz, tx, ty, tz.
    """

    mass: NDArray[np.float64]
    stiffness: NDArray[np.float64]
    dofs: tuple[Dof, ...]


def assemble(structure: Structure, *, progress: Progress | None = None) -> SystemMatrices:
    """Assemble mass and stiffness matrices; the anchor's rows and columns are dropped.

    Stage assemble logs through the ``connflex_core`` logger and reports to ``progress``.
    """
    with stage("assemble", progress) as current:
        system = _assemble(structure)
        current.count(dofs=len(system.dofs), beams=len(structure.beams))
    return system


def applied_loads(
    system: SystemMatrices, loads: Mapping[int, Sequence[float]] | None
) -> NDArray[np.float64]:
    """Generalized forces on the rows of ``system`` from ``loads``, a mapping of body index
    to six forces and moments at its center of mass; None or empty applies nothing."""
    rows = {body: row for row, (body, name) in enumerate(system.dofs) if name == "ux"}
    force = np.zeros(len(system.dofs))
    for body, load in (loads or {}).items():
        if body not in rows:
            raise ValueError(f"body {body} is the anchor or does not exist")
        if len(load) != 6:
            raise ValueError(f"load on body {body} must have six components")
        force[rows[body] : rows[body] + 6] = load
    return force


def anchored(structure: Structure, system: SystemMatrices) -> NDArray[np.bool_]:
    """Which rows of ``system`` belong to bodies that a chain of beams joins to the
    anchor; the others cannot balance a force. None at all is an error."""
    joined = nx.node_connected_component(structure.graph, structure.anchor.index)
    kept = np.array([body in joined for body, _ in system.dofs], dtype=bool)
    if not kept.any():
        raise ValueError("no movable body is joined to the anchor by beams")
    return kept


def _assemble(structure: Structure) -> SystemMatrices:
    movable = [body for body in structure.bodies if not body.anchor]
    offset = {body.index: 6 * position for position, body in enumerate(movable)}
    size = 6 * len(movable)
    mass = np.zeros((size, size))
    stiffness = np.zeros((size, size))
    for body in movable:
        row = offset[body.index]
        mass[row : row + 6, row : row + 6] = body.inertial.matrix
    centers = {body.index: body.inertial.center for body in structure.bodies}
    material = structure.stack.device.material
    for beam in structure.beams:
        a, b = beam.bodies
        k = effective_stiffness(beam, material, centers[a], centers[b])
        for i, body_i in enumerate((a, b)):
            for j, body_j in enumerate((a, b)):
                if body_i in offset and body_j in offset:
                    ri, rj = offset[body_i], offset[body_j]
                    stiffness[ri : ri + 6, rj : rj + 6] += k[6 * i : 6 * i + 6, 6 * j : 6 * j + 6]
    dofs = tuple((body.index, name) for body in movable for name in DOF_NAMES)
    return SystemMatrices(mass=mass, stiffness=stiffness, dofs=dofs)

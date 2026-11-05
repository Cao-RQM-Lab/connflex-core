"""Static equilibrium under electrode potentials."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

import networkx as nx
import numpy as np
import scipy.linalg
from numpy.typing import NDArray

from connflex_core.mechanics.assembly import Dof, assemble
from connflex_core.params import StaticOptions
from connflex_core.progress import Progress, stage

if TYPE_CHECKING:
    from connflex_core.electrostatics.capacitance import CapacitanceModel
    from connflex_core.model import Structure

Array = NDArray[np.float64]

HALVINGS = 12


@dataclass(frozen=True, eq=False)
class Equilibrium:
    """Where a structure rests under electrode potentials, equations (115) to (118).

    Attributes:
        displacement: Displacement of every degree of freedom in ``dofs`` from rest, in
            meters and radians.
        dofs: (body index, degree of freedom) of each entry, as for the assembled
            mechanical matrices.
        potentials: Potential of every electrode, in volts.
        charge: Charge on every electrode, in coulombs.
        force: Electrostatic generalized forces at ``displacement``.
        stiffness: Tangent stiffness there: mechanical plus electrostatic.
        softest: Smallest eigenvalue of ``stiffness`` against the mass matrix, the
            squared angular frequency of the softest motion under these potentials.
        stable: Whether ``softest`` is positive.
        iterations: Newton iterations taken.
        loose: Bodies that no chain of beams joins to the anchor. Nothing can balance a
            force on them, so they are kept at rest and left out of ``softest``.
        model: Capacitance model at ``displacement``.
    """

    displacement: Array
    dofs: tuple[Dof, ...]
    potentials: Array
    charge: Array
    force: Array
    stiffness: Array
    softest: float
    stable: bool
    iterations: int
    loose: tuple[int, ...]
    model: CapacitanceModel


def equilibrium(
    structure: Structure,
    model: CapacitanceModel,
    potentials: Mapping[int, float],
    *,
    loads: Mapping[int, Sequence[float]] | None = None,
    options: StaticOptions | None = None,
    progress: Progress | None = None,
) -> Equilibrium:
    """Displacement at which beam forces balance electrostatic forces.

    ``model`` is the capacitance model of ``structure`` at rest. ``potentials`` maps a
    net to its potential in volts; a net it does not name is at zero. ``loads`` maps a
    body index to six generalized forces at its center of mass, as for ``static``.
    Newton iteration starts at rest, equation (116), and evaluates the capacitance model
    at every iterate. No equilibrium within ``max_iterations`` is an error: potentials
    may be beyond pull-in, or the plates method may have a kink there. Stages assemble
    and equilibrium log through the
    ``connflex_core`` logger and report to ``progress``.
    """
    options = options or StaticOptions()
    if model.rest.structure is not structure:
        raise ValueError("capacitance model was built from another structure")
    volts = np.zeros(structure.nets.count)
    for net, value in potentials.items():
        if not 0 <= net < len(volts):
            raise ValueError(f"net {net} does not exist")
        volts[net] = value
    with stage("assemble", progress) as current:
        system = assemble(structure)
        current.count(dofs=len(system.dofs), beams=len(structure.beams))
    joined = nx.node_connected_component(structure.graph, structure.anchor.index)
    kept = np.array([body in joined for body, _ in system.dofs], dtype=bool)
    if not kept.any():
        raise ValueError("no movable body is joined to the anchor by beams")
    rows = {body: row for row, (body, name) in enumerate(system.dofs) if name == "ux"}
    applied = np.zeros(len(system.dofs))
    for body, load in (loads or {}).items():
        if body not in rows:
            raise ValueError(f"body {body} is the anchor or does not exist")
        if len(load) != 6:
            raise ValueError(f"load on body {body} must have six components")
        applied[rows[body] : rows[body] + 6] = load

    held = volts[model.net]
    with stage("equilibrium", progress) as current:
        found = _settle(model, system.stiffness, kept, held, applied, options)
        if found is None:
            raise ValueError(
                f"no equilibrium in {options.max_iterations} iterations: potentials may be "
                f'beyond pull-in, or forces of sidewall method "plates" may jump where '
                f"corners of two bodies line up"
            )
        here, iterations = found
        current.count(iterations=iterations)
        tangent = system.stiffness + here.stiffness(held)
        both = np.ix_(kept, kept)
        softest = float(scipy.linalg.eigh(tangent[both], system.mass[both], eigvals_only=True)[0])
    return Equilibrium(
        displacement=here.displacement,
        dofs=system.dofs,
        potentials=held,
        charge=here.matrix @ held,
        force=here.force(held),
        stiffness=tangent,
        softest=softest,
        stable=softest > 0,
        iterations=iterations,
        loose=tuple(sorted({body for body, _ in system.dofs if body not in joined})),
        model=here,
    )


def _settle(
    model: CapacitanceModel,
    spring: Array,
    kept: NDArray[np.bool_],
    held: Array,
    applied: Array,
    options: StaticOptions,
) -> tuple[CapacitanceModel, int] | None:
    """Newton iteration from the displacement of ``model``; None when it does not converge.

    Only degrees of freedom marked in ``kept`` move.
    """
    here = model
    both = np.ix_(kept, kept)
    settling = False
    for iteration in range(1, options.max_iterations + 1):
        place = here.displacement
        unbalance = spring @ place - here.force(held) - applied
        tangent = spring + here.stiffness(held)
        step = np.zeros(len(place))
        try:
            step[kept] = scipy.linalg.solve(tangent[both], -unbalance[kept], assume_a="sym")
        except scipy.linalg.LinAlgError:
            return None
        whole = True
        for _ in range(HALVINGS):
            try:
                there = model.at(place + step)
                break
            except ValueError:
                step, whole = step / 2, False
        else:
            return None
        here = there
        size = float(step @ spring @ step)
        strain = float(here.displacement @ spring @ here.displacement)
        if whole and (settling or size <= options.tolerance**2 * strain):
            return here, iteration
        settling = whole and size <= options.tolerance * strain
    return None

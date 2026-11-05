"""Static equilibrium under electrode potentials."""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
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
PATIENCE = 4


@dataclass(frozen=True, eq=False)
class Equilibrium:
    """Where a structure rests under electrode potentials, equations (115) to (123).

    Attributes:
        scale: Fraction of the potentials asked for at which this equilibrium holds: 1
            when they were reached, else the largest fraction with a stable equilibrium.
        pull_in: Fraction of the potentials asked for at which no stable equilibrium was
            found any more, within ``pull_in_resolution`` above ``scale``. None when the
            potentials were reached.
        limit: Why potentials stopped rising: ``"contact"`` when the equilibrium at
            ``pull_in`` closes a gap to ``contact_fraction`` of its width,
            ``"unstable"`` when it does not hold against a push, ``"no equilibrium"``
            when none was found. None when the potentials were reached.
        displacement: Displacement of every degree of freedom in ``dofs`` from rest, in
            meters and radians.
        dofs: (body index, degree of freedom) of each entry, as for the assembled
            mechanical matrices.
        potentials: Potential of every electrode at ``scale``, in volts, those of
            floating nets as their charge sets them.
        charge: Charge on every electrode, in coulombs.
        force: Electrostatic generalized forces at ``displacement``.
        stiffness: Tangent stiffness there: mechanical plus electrostatic, with the
            part that charge held on floating nets adds, equation (123).
        softest: Smallest eigenvalue of ``stiffness`` against the mass matrix, the
            squared angular frequency of the softest motion under these potentials.
        stable: Whether ``softest`` is positive.
        iterations: Newton iterations taken, over all steps.
        loose: Bodies that no chain of beams joins to the anchor. Nothing can balance a
            force on them, so they are kept at rest and left out of ``softest``.
        model: Capacitance model at ``displacement``.
        ramp_scales: Fractions of the potentials at which a stable equilibrium was found
            on the way, ascending from zero to ``scale``.
        ramp_displacements: Displacement at each of them, one row each.
        ramp_softest: ``softest`` at each of them.
        ramp_clearance: ``CapacitanceModel.clearance`` at each of them.
    """

    scale: float
    pull_in: float | None
    limit: str | None
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
    ramp_scales: Array
    ramp_displacements: Array
    ramp_softest: Array
    ramp_clearance: Array


def equilibrium(
    structure: Structure,
    model: CapacitanceModel,
    potentials: Mapping[int, float],
    *,
    charges: Mapping[int, float] | None = None,
    loads: Mapping[int, Sequence[float]] | None = None,
    options: StaticOptions | None = None,
    progress: Progress | None = None,
) -> Equilibrium:
    """Displacement at which beam forces balance electrostatic forces.

    ``model`` is the capacitance model of ``structure`` at rest. ``potentials`` maps a
    net to its potential in volts. ``charges`` maps a floating net to the charge it
    holds, in coulombs; its potential follows from equation (122). A net in neither is
    at zero volts. ``loads`` maps a body index to six generalized forces at its center
    of mass, as for ``static``. Potentials and charges rise together from zero.
    Potentials rise from zero in ``steps`` and the equilibrium is followed by Newton
    iteration, equation (116), which evaluates the capacitance model at every iterate.
    Where no stable equilibrium is found the step is halved, and the result is the last
    stable equilibrium with the fraction of the potentials at which stability was lost,
    equation (119). Stages assemble and equilibrium log through the ``connflex_core``
    logger and report to ``progress``; equilibrium advances once per step.
    """
    options = options or StaticOptions()
    if model.rest.structure is not structure:
        raise ValueError("capacitance model was built from another structure")
    volts = np.zeros(structure.nets.count)
    stored = np.zeros(structure.nets.count)
    for values, given in ((volts, potentials), (stored, charges or {})):
        for net, value in given.items():
            if not 0 <= net < len(volts):
                raise ValueError(f"net {net} does not exist")
            values[net] = value
    floating = np.array(sorted(charges or {}), dtype=np.intp)
    if set(floating.tolist()) & set(potentials):
        raise ValueError("a net cannot have both a potential and a charge")
    driven = np.setdiff1d(np.arange(len(volts)), floating)
    tie = (model.net[:, None] == np.arange(len(volts))).astype(float)
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

    both = np.ix_(kept, kept)

    def drive(scale: float) -> Callable[[CapacitanceModel], tuple[Array, Array]]:
        """Electrode potentials at a displacement, and stiffness that held charge adds."""

        def at(here: CapacitanceModel) -> tuple[Array, Array]:
            level = scale * volts
            added = np.zeros_like(system.stiffness)
            if floating.size:
                between = tie.T @ here.matrix @ tie
                within = between[np.ix_(floating, floating)]
                try:
                    factor = scipy.linalg.cho_factor(within)
                except scipy.linalg.LinAlgError as exc:
                    raise ValueError(
                        "a floating net has no capacitance to a net with a potential"
                    ) from exc
                across = between[np.ix_(floating, driven)] @ level[driven]
                level[floating] = scipy.linalg.cho_solve(factor, scale * stored[floating] - across)
                pull = (tie.T @ np.einsum("abk,b->ak", here.gradient, level[here.net]))[floating]
                added = pull.T @ scipy.linalg.cho_solve(factor, pull)
            return level[here.net], added

        return at

    def softness(here: CapacitanceModel, scale: float) -> float:
        held, added = drive(scale)(here)
        tangent = system.stiffness + here.stiffness(held) + added
        return float(scipy.linalg.eigh(tangent[both], system.mass[both], eigvals_only=True)[0])

    start, iterations = _settle(model, system.stiffness, kept, drive(0.0), applied, options)
    if start is None:
        raise ValueError("no equilibrium under mechanical loads alone")
    here = start
    scale, pull_in, limit = 0.0, None, None
    ramp = [(scale, here.displacement, softness(here, scale), here.clearance)]
    from_here = False
    with stage("equilibrium", progress, total=options.steps) as current:
        while scale < 1:
            if pull_in is None:
                target = min((np.floor(scale * options.steps + 1e-9) + 1) / options.steps, 1.0)
            elif pull_in - scale > options.pull_in_resolution:
                target = (scale + pull_in) / 2
            elif from_here:
                break
            else:
                target = pull_in
            found, taken = _settle(here, system.stiffness, kept, drive(target), applied, options)
            softest = softness(found, target) if found else 0.0
            iterations += taken
            if found is not None and softest > 0 and found.clearance > options.contact_fraction:
                if pull_in is None:
                    current.advance()
                elif target == pull_in:
                    pull_in = limit = None
                here, scale, from_here = found, target, False
                ramp.append((scale, here.displacement, softest, here.clearance))
            else:
                pull_in, from_here = target, True
                limit = (
                    "no equilibrium" if found is None else "unstable" if softest <= 0 else "contact"
                )
        current.count(iterations=iterations, steps=len(ramp) - 1)
    reached, added = drive(scale)(here)
    return Equilibrium(
        scale=scale,
        pull_in=pull_in,
        limit=limit,
        displacement=here.displacement,
        dofs=system.dofs,
        potentials=reached,
        charge=here.matrix @ reached,
        force=here.force(reached),
        stiffness=system.stiffness + here.stiffness(reached) + added,
        softest=ramp[-1][2],
        stable=ramp[-1][2] > 0,
        iterations=iterations,
        loose=tuple(sorted({body for body, _ in system.dofs if body not in joined})),
        model=here,
        ramp_scales=np.array([entry[0] for entry in ramp]),
        ramp_displacements=np.array([entry[1] for entry in ramp]),
        ramp_softest=np.array([entry[2] for entry in ramp]),
        ramp_clearance=np.array([entry[3] for entry in ramp]),
    )


def _settle(
    model: CapacitanceModel,
    spring: Array,
    kept: NDArray[np.bool_],
    drive: Callable[[CapacitanceModel], tuple[Array, Array]],
    applied: Array,
    options: StaticOptions,
) -> tuple[CapacitanceModel | None, int]:
    """Newton iteration from the displacement of ``model``, with the iterations it took.

    ``drive`` gives electrode potentials at a displacement and the stiffness that held
    charge adds. Only degrees of freedom marked in ``kept`` move. The model is None when iteration
    does not converge; it is given up after ``PATIENCE`` steps in a row that are halved
    or are no smaller than the smallest so far.
    """
    here = model
    both = np.ix_(kept, kept)
    settling = False
    smallest, idle = np.inf, 0
    for iteration in range(1, options.max_iterations + 1):
        place = here.displacement
        held, added = drive(here)
        unbalance = spring @ place - here.force(held) - applied
        tangent = spring + here.stiffness(held) + added
        step = np.zeros(len(place))
        try:
            step[kept] = scipy.linalg.solve(tangent[both], -unbalance[kept], assume_a="sym")
        except scipy.linalg.LinAlgError:
            return None, iteration
        whole = True
        for _ in range(HALVINGS):
            try:
                there = model.at(place + step)
                break
            except ValueError:
                step, whole = step / 2, False
        else:
            return None, iteration
        here = there
        size = float(step @ spring @ step)
        strain = float(here.displacement @ spring @ here.displacement)
        if whole and (settling or size <= options.tolerance**2 * strain):
            return here, iteration
        settling = whole and size <= options.tolerance * strain
        idle = 0 if whole and size < smallest else idle + 1
        smallest = min(smallest, size)
        if idle == PATIENCE:
            return None, iteration
    return None, options.max_iterations

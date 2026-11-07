"""Static equilibrium under electrode potentials."""

from __future__ import annotations

import warnings
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

import networkx as nx
import numpy as np
import scipy.linalg
from numpy.typing import ArrayLike, NDArray

from connflex_core.mechanics.assembly import Dof, assemble
from connflex_core.params import StaticOptions
from connflex_core.progress import Progress, stage

if TYPE_CHECKING:
    from connflex_core.electrostatics.capacitance import CapacitanceModel
    from connflex_core.model import Structure

Array = NDArray[np.float64]

HALVINGS = 12


class NoEquilibriumError(ValueError):
    """Newton iteration found no equilibrium from the start it was given."""


@dataclass(frozen=True, eq=False)
class Equilibrium:
    """Where a structure rests under electrode potentials, equations (120) to (128).

    Attributes:
        displacement: Displacement of every degree of freedom in ``dofs`` from rest, in
            meters and radians.
        dofs: (body index, degree of freedom) of each entry, as for the assembled
            mechanical matrices.
        potentials: Potential of every electrode, in volts, those of floating nets as
            their charge sets them.
        charge: Charge on every electrode, in coulombs.
        force: Electrostatic generalized forces at ``displacement``.
        stiffness: Tangent stiffness there: mechanical plus electrostatic, with the
            part that charge held on floating nets adds, equation (128).
        softest: Smallest eigenvalue of ``stiffness`` against the mass matrix, the
            squared angular frequency of the softest motion under these potentials.
        margin: Smallest eigenvalue of ``stiffness`` against the mechanical stiffness:
            the fraction of beam stiffness left in the direction the field softens
            most, equation (124). 1 without potentials, 0 at pull-in, negative when the
            equilibrium does not hold against a push.
        stable: Whether ``softest`` is positive.
        iterations: Newton iterations taken.
        loose: Bodies that no chain of beams joins to the anchor. Nothing can balance a
            force on them, so they are kept at rest and left out of ``softest``.
        floating: Nets that held a fixed charge, ascending.
        model: Capacitance model at ``displacement``, whose ``clearance`` tells how far
            the narrowest gap between two bodies has closed.
    """

    displacement: Array
    dofs: tuple[Dof, ...]
    potentials: Array
    charge: Array
    force: Array
    stiffness: Array
    softest: float
    margin: float
    stable: bool
    iterations: int
    loose: tuple[int, ...]
    floating: tuple[int, ...]
    model: CapacitanceModel


def equilibrium(
    structure: Structure,
    model: CapacitanceModel,
    potentials: Mapping[int, float],
    *,
    charges: Mapping[int, float] | None = None,
    loads: Mapping[int, Sequence[float]] | None = None,
    start: Equilibrium | ArrayLike | None = None,
    options: StaticOptions | None = None,
    workers: int | None = None,
    progress: Progress | None = None,
) -> Equilibrium:
    """Displacement at which beam forces balance electrostatic forces, one operating point.

    ``model`` is the capacitance model of ``structure`` at rest. ``potentials`` maps a
    net to its potential in volts. ``charges`` maps a floating net to the charge it
    holds, in coulombs; its potential follows from equation (127). A net in neither is
    at zero volts. ``loads`` maps a body index to six generalized forces at its center
    of mass, as for ``static``. Newton iteration, equation (121), starts at rest, or at
    ``start``: an earlier ``Equilibrium`` of this structure or a displacement of every
    degree of freedom. It finds the equilibrium near its start, so which of several it
    finds depends on where it starts; ``stable`` and ``margin`` tell whether that one
    holds. Iteration that has not converged within ``max_iterations``, whose steps stop
    shrinking, whose tangent is singular, or that meets no finite force (theory,
    Implementation choices) raises ``NoEquilibriumError``.
    ``workers`` threads evaluate the model, as many as at rest by default. Stages
    assemble and equilibrium log through the ``connflex_core`` logger and report to
    ``progress``.
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
    origin = _origin(model, start, kept, workers)

    def drive(here: CapacitanceModel) -> tuple[Array, Array]:
        """Electrode potentials at a displacement, and stiffness that held charge adds."""
        level = volts.copy()
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
            level[floating] = scipy.linalg.cho_solve(factor, stored[floating] - across)
            pull = (tie.T @ np.einsum("abk,b->ak", here.gradient, level[here.net]))[floating]
            added = pull.T @ scipy.linalg.cho_solve(factor, pull)
        return level[here.net], added

    with stage("equilibrium", progress) as current:
        here, iterations = _settle(origin, system.stiffness, kept, drive, applied, options, workers)
        current.count(iterations=iterations)
        if here is None:
            raise NoEquilibriumError(
                f"no equilibrium found in {iterations} iterations from this start: the "
                f"potentials may be beyond pull-in, or the start too far from an "
                f"equilibrium; start from an equilibrium at nearer potentials"
            )
    if model.rest.params.surrogate_tolerance is not None:
        here = here.at(here.displacement, exact=True, workers=workers)
    reached, added = drive(here)
    both = np.ix_(kept, kept)
    tangent = system.stiffness + here.stiffness(reached) + added
    softest = float(scipy.linalg.eigh(tangent[both], system.mass[both], eigvals_only=True)[0])
    margin = float(scipy.linalg.eigh(tangent[both], system.stiffness[both], eigvals_only=True)[0])
    return Equilibrium(
        displacement=here.displacement,
        dofs=system.dofs,
        potentials=reached,
        charge=here.matrix @ reached,
        force=here.force(reached),
        stiffness=tangent,
        softest=softest,
        margin=margin,
        stable=softest > 0,
        iterations=iterations,
        loose=tuple(sorted({body for body, _ in system.dofs if body not in joined})),
        floating=tuple(floating.tolist()),
        model=here,
    )


def _origin(
    model: CapacitanceModel,
    start: Equilibrium | ArrayLike | None,
    kept: NDArray[np.bool_],
    workers: int | None,
) -> CapacitanceModel:
    """Capacitance model where Newton iteration starts."""
    if start is None:
        return model
    if isinstance(start, Equilibrium):
        rest, other = model.rest, start.model.rest
        if other.structure is not rest.structure or other.params != rest.params:
            raise ValueError("start is an equilibrium of another structure or model")
        return start.model
    place = np.asarray(start, dtype=float)
    if place.shape != (len(model.dofs),):
        raise ValueError(
            f"start must hold one value for each of {len(model.dofs)} degrees of freedom"
        )
    return model.at(np.where(kept, place, 0.0), workers=workers)


def _settle(
    model: CapacitanceModel,
    spring: Array,
    kept: NDArray[np.bool_],
    drive: Callable[[CapacitanceModel], tuple[Array, Array]],
    applied: Array,
    options: StaticOptions,
    workers: int | None,
) -> tuple[CapacitanceModel | None, int]:
    """Newton iteration from the displacement of ``model``, with the iterations it took.

    ``drive`` gives electrode potentials at a displacement and the stiffness that held
    charge adds. Only degrees of freedom marked in ``kept`` move. The model is None when
    iteration does not converge: after ``patience`` steps in a row that are halved or are
    no smaller than the smallest so far, after ``HALVINGS`` halvings of one step, on a
    singular tangent, or on a force or stiffness that is not finite.
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
        if not (np.isfinite(unbalance[kept]).all() and np.isfinite(tangent[both]).all()):
            return None, iteration
        step = np.zeros(len(place))
        try:
            with warnings.catch_warnings():  # singular along a direction without unbalance
                warnings.simplefilter("ignore", scipy.linalg.LinAlgWarning)
                step[kept] = scipy.linalg.solve(tangent[both], -unbalance[kept], assume_a="sym")
        except scipy.linalg.LinAlgError:
            return None, iteration
        whole = True
        for _ in range(HALVINGS):
            try:
                there = model.at(place + step, workers=workers)
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
        if idle == options.patience:
            return None, iteration
    return None, options.max_iterations

"""Motion, potentials, and charges over time under sources."""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from itertools import pairwise
from typing import TYPE_CHECKING, Any

import networkx as nx
import numpy as np
import scipy.integrate
import scipy.linalg
import scipy.optimize
from numpy.typing import ArrayLike, NDArray

from connflex_core.coupled.static import Equilibrium
from connflex_core.electrical.drivers import Source
from connflex_core.electrical.network import electrode_indices, network_nodes
from connflex_core.mechanics.assembly import Dof, assemble
from connflex_core.params import TransientOptions
from connflex_core.progress import Progress, stage

if TYPE_CHECKING:
    from connflex_core.electrostatics.capacitance import CapacitanceModel
    from connflex_core.model import Structure

Array = NDArray[np.float64]

RETRIES = 8
CACHE = 16


class _ClosedError(Exception):
    """A trial step brought two bodies into contact."""


@dataclass(frozen=True)
class Energy:
    """Energy ledger of a transient at each output time, in joules, equation (145).

    Attributes:
        kinetic: Kinetic energy of the bodies.
        strain: Strain energy in the beams.
        electrostatic: Energy in the field between electrodes.
        damped: Energy that damping has taken since the start.
        resistive: Energy that resistance of links and sources has taken since the start.
        supplied: Energy that sources and loads have put in since the start.
    """

    kinetic: Array
    strain: Array
    electrostatic: Array
    damped: Array
    resistive: Array
    supplied: Array

    @property
    def imbalance(self) -> Array:
        """Stored energy gained plus energy taken, less energy supplied: zero but for
        the error of integration."""
        stored = self.kinetic + self.strain + self.electrostatic
        imbalance: Array = stored - stored[0] + self.damped + self.resistive - self.supplied
        return imbalance


@dataclass(frozen=True, eq=False)
class Transient:
    """Response of a structure to sources over time, equations (139) to (147).

    Attributes:
        times: Output times in seconds: those asked for, up to the end of the run, and
            the time of a contact event last when one ended it.
        displacement: Displacement of every degree of freedom in ``dofs`` at each time.
        velocity: Velocity of every degree of freedom at each time.
        potentials: Potential of every electrode at each time.
        charge: Charge on every electrode at each time.
        current: Current that each source drives into its electrode at each time, in
            the order of the sources.
        clearance: Clearance of equation (126) at each time.
        energy: Energy ledger at each time.
        contact: Electrodes across the gap whose closing to ``contact_clearance`` ended
            the run. None when the run reached its last time.
        dofs: (body index, degree of freedom) of each column of ``displacement``.
        loose: Bodies that no chain of beams joins to the anchor, kept at rest.
        floating: Nets that kept a fixed charge.
        steps: Time steps the integrator accepted.
        evaluations: Evaluations of the capacitance model at displaced bodies.
        model: Capacitance model at the last time.
    """

    times: Array
    displacement: Array
    velocity: Array
    potentials: Array
    charge: Array
    current: Array
    clearance: Array
    energy: Energy
    contact: tuple[int, int] | None
    dofs: tuple[Dof, ...]
    loose: tuple[int, ...]
    floating: tuple[int, ...]
    steps: int
    evaluations: int
    model: CapacitanceModel


@dataclass
class _State:
    """Everything one evaluation at a time and state gives."""

    here: CapacitanceModel
    potentials: Array
    electrodes: Array
    factor: Any
    joined: Array
    rates: Array
    current: Array


@dataclass(eq=False)
class _Run:
    """Network, mechanics, and sources of a run, with its right side and Jacobian."""

    model: CapacitanceModel
    sources: Sequence[Source]
    tie: Array
    ideal: NDArray[np.intp]
    dynamic: NDArray[np.intp]
    unknown: NDArray[np.intp]
    floating: Array
    conductance: Array
    feed: Array
    kept: NDArray[np.bool_]
    mass: Array
    spring: Array
    damping: Array
    applied: Array
    duration: float
    evaluations: int = 0
    cache: dict[bytes, CapacitanceModel] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.factor = scipy.linalg.cho_factor(self.mass)
        self.count = int(self.kept.sum())
        self.dofs = len(self.kept)
        self.held = np.array([source.resistance == 0 for source in self.sources], dtype=bool)
        self.slowed = scipy.linalg.cho_solve(self.factor, self.damping)
        self.total = self.conductance + np.diag(self.feed)

    def at(self, place: Array) -> CapacitanceModel:
        key = place.tobytes()
        if key not in self.cache:
            full = np.zeros(self.dofs)
            full[self.kept] = place
            try:
                self.cache[key] = self.model.at(full)
            except ValueError as exc:
                raise _ClosedError from exc
            self.evaluations += 1
            if len(self.cache) > CACHE:
                del self.cache[next(iter(self.cache))]
        return self.cache[key]

    def emf(self, time: float) -> Array:
        return np.array([_value(source.waveform, time) for source in self.sources])

    def slope(self, time: float) -> Array:
        return np.array([_slope(source.waveform, time, self.duration) for source in self.sources])

    def state(self, time: float, y: Array) -> _State:
        """Potentials, rates of every state, and source currents, equations (139) to (144)."""
        n, d = self.count, len(self.dynamic)
        place, speed, held = y[:n], y[n : 2 * n], y[2 * n : 2 * n + d]
        here = self.at(place)
        joined = self.tie.T @ here.matrix @ self.tie
        emf = self.emf(time)
        nodes = np.zeros(len(joined))
        nodes[self.ideal] = emf[self.held]
        unknown, ideal = self.unknown, self.ideal
        factor = None
        if unknown.size:
            factor = scipy.linalg.cho_factor(joined[np.ix_(unknown, unknown)])
            given = np.concatenate([held, self.floating]) - joined[np.ix_(unknown, ideal)] @ (
                nodes[ideal]
            )
            nodes[unknown] = scipy.linalg.cho_solve(factor, given)
        electrodes = self.tie @ nodes
        force = here.force(electrodes)[self.kept]
        acceleration = scipy.linalg.cho_solve(
            self.factor, force + self.applied - self.spring @ place - self.damping @ speed
        )
        driven = len(self.conductance)
        source_emf = np.zeros(driven)
        source_emf[: len(self.sources)] = emf
        flowing = -self.total @ nodes[:driven] + self.feed * source_emf
        charging = flowing[self.dynamic]

        velocity = np.zeros(self.dofs)
        velocity[self.kept] = speed
        changing = self.tie.T @ np.einsum("abk,k->ab", here.gradient, velocity) @ self.tie
        drift = changing @ nodes
        rising = np.zeros(len(joined))
        rising[ideal] = self.slope(time)[self.held]
        if factor is not None:
            pushed = np.concatenate([charging, np.zeros(len(self.floating))])
            pushed -= drift[unknown] + joined[np.ix_(unknown, ideal)] @ rising[ideal]
            rising[unknown] = scipy.linalg.cho_solve(factor, pushed)
        held_current = (
            drift[ideal] + joined[ideal] @ rising + (self.conductance @ nodes[:driven])[ideal]
        )
        current = np.zeros(len(self.sources))
        current[self.held] = held_current
        resistive = ~self.held
        current[resistive] = (self.feed * (source_emf - nodes[:driven]))[: len(self.sources)][
            resistive
        ]

        damped = float(speed @ self.damping @ speed)
        across = nodes[:driven]
        lost = float(across @ self.conductance @ across) + float(
            self.feed @ (source_emf - across) ** 2
        )
        supplied = float(emf @ current) + float(self.applied @ speed)
        rates = np.concatenate([speed, acceleration, charging, [damped, lost, supplied]])
        return _State(here, nodes, electrodes, factor, joined, rates, current)

    def fun(self, time: float, y: Array) -> Array:
        return self.state(time, y).rates

    def jac(self, time: float, y: Array) -> Array:
        """Jacobian of the rates in the states, equation (143)."""
        now = self.state(time, y)
        n, d = self.count, len(self.dynamic)
        here, unknown = now.here, self.unknown
        stiffness = here.stiffness(now.electrodes)[np.ix_(self.kept, self.kept)]
        pull = (self.tie.T @ np.einsum("abk,b->ak", here.gradient, now.electrodes))[:, self.kept]
        size = 2 * n + d + 3
        matrix = np.zeros((size, size))
        matrix[:n, n : 2 * n] = np.eye(n)
        matrix[n : 2 * n, n : 2 * n] = -self.slowed
        pulled = -stiffness - self.spring
        if unknown.size:
            reach = pull[unknown]
            moved = scipy.linalg.cho_solve(now.factor, reach)
            charged = scipy.linalg.cho_solve(now.factor, np.eye(len(unknown))[:, :d])
            pulled -= reach.T @ moved
            matrix[n : 2 * n, 2 * n : 2 * n + d] = scipy.linalg.cho_solve(
                self.factor, reach.T @ charged
            )
            leak = self.total[np.ix_(self.dynamic, self.dynamic)]
            matrix[2 * n : 2 * n + d, :n] = leak @ moved[:d]
            matrix[2 * n : 2 * n + d, 2 * n : 2 * n + d] = -leak @ charged[:d]
        matrix[n : 2 * n, :n] = scipy.linalg.cho_solve(self.factor, pulled)
        return matrix


def transient(
    structure: Structure,
    model: CapacitanceModel,
    sources: Sequence[Source],
    times: ArrayLike,
    options: TransientOptions,
    *,
    charges: Mapping[int, float] | None = None,
    loads: Mapping[int, Sequence[float]] | None = None,
    start: Equilibrium | ArrayLike | None = None,
    progress: Progress | None = None,
) -> Transient:
    """Motion, potentials, and charges of ``structure`` over ``times`` under ``sources``.

    ``model`` is the capacitance model of ``structure`` at rest. Each source drives one
    electrode; electrodes of its net that links with resistance set apart from it are
    reached through them. ``charges`` maps a floating net to the charge it keeps. A net
    with neither is held at zero volts on every electrode. ``loads`` are constant forces
    on bodies, as for ``static``. The run starts at ``times[0]`` from rest or from the
    displacement of ``start``, at zero velocity, with the potentials that the sources
    hold without current at that time, and ends at ``times[-1]`` or at contact. Stages
    assemble and transient log through the ``connflex_core`` logger and report to
    ``progress``, which advances once per output time.
    """
    if model.rest.structure is not structure:
        raise ValueError("capacitance model was built from another structure")
    stamps = np.atleast_1d(np.asarray(times, dtype=float))
    if stamps.ndim != 1 or len(stamps) < 2 or np.any(np.diff(stamps) <= 0):
        raise ValueError("times must hold at least two strictly ascending values")
    if not np.all(np.isfinite(stamps)):
        raise ValueError("times must be finite")
    stored = dict(charges or {})
    for net in stored:
        if not 0 <= net < structure.nets.count:
            raise ValueError(f"net {net} does not exist")
    listed: Sequence[object] = sources
    for source in listed:
        if not isinstance(source, Source):
            raise TypeError(f"a source is a Source, not {type(source).__name__} {source!r}")
    places = electrode_indices(
        structure, [source.electrode for source in sources], stored, "source"
    )
    with stage("assemble", progress) as current:
        system = assemble(structure)
        current.count(dofs=len(system.dofs), beams=len(structure.beams))
    joined_to_anchor = nx.node_connected_component(structure.graph, structure.anchor.index)
    kept = np.array([body in joined_to_anchor for body, _ in system.dofs], dtype=bool)
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
    origin = _origin(model, start, kept)

    nodes, hidden, conductance = network_nodes(structure, stored, places, "source")
    tie = np.zeros((len(structure.electrodes), len(nodes) + len(hidden)))
    for column, node in enumerate([*nodes, *hidden]):
        tie[list(node), column] = 1.0
    held = np.array([source.resistance == 0 for source in sources], dtype=bool)
    ideal = np.flatnonzero(held).astype(np.intp)
    dynamic = np.setdiff1d(np.arange(len(nodes)), ideal).astype(np.intp)
    unknown = np.concatenate([dynamic, np.arange(len(nodes), len(nodes) + len(hidden))])
    feed = np.zeros(len(nodes))
    for k, source in enumerate(sources):
        if source.resistance:
            feed[k] = 1 / source.resistance
    floating = np.array(
        [stored[structure.electrodes[node[0]].net] for node in hidden], dtype=np.float64
    )

    both = np.ix_(kept, kept)
    mass, spring = system.mass[both], system.stiffness[both]
    squared, shapes = scipy.linalg.eigh(spring, mass)
    natural = np.sqrt(squared)
    weighted = mass @ shapes
    damping = weighted @ np.diag(natural / options.quality_factor) @ weighted.T
    run = _Run(
        model=model,
        sources=sources,
        tie=tie,
        ideal=ideal,
        dynamic=dynamic,
        unknown=unknown.astype(np.intp),
        floating=floating,
        conductance=conductance,
        feed=feed,
        kept=kept,
        mass=mass,
        spring=spring,
        damping=damping,
        applied=applied[kept],
        duration=float(stamps[-1] - stamps[0]),
    )
    start_place = origin[kept]
    y = np.concatenate(
        [start_place, np.zeros(run.count), _rest_charge(run, start_place, stamps[0]), np.zeros(3)]
    )
    scale = _scales(structure, run, system.dofs, natural, stamps, options.tolerance)
    breaks = sorted(
        {
            float(point)
            for source in sources
            for point in getattr(source.waveform, "breakpoints", ())
            if stamps[0] < point < stamps[-1]
        }
    )
    frames: list[tuple[float, Array]] = []
    contact: tuple[int, int] | None = None
    steps = 0
    with stage("transient", progress, total=len(stamps)) as current:
        frames.append((float(stamps[0]), y))
        current.advance()
        pending = 1
        edges = [float(stamps[0]), *breaks, float(stamps[-1])]
        for left, right in pairwise(edges):
            time = left
            limit = options.max_step or np.inf
            solver = _solver(run, time, y, right, options.tolerance, scale, limit)
            tries = 0
            while solver.status == "running":
                before, state_before = solver.t, solver.y.copy()
                last = solver.step_size or (right - left)
                try:
                    solver.step()
                except _ClosedError:
                    tries += 1
                    if tries > RETRIES:
                        contact = run.at(state_before[: run.count]).closest or (-1, -1)
                        frames.append((before, state_before))
                        break
                    limit = last / 4**tries
                    solver = _solver(
                        run, before, state_before, right, options.tolerance, scale, limit
                    )
                    continue
                tries = 0
                if solver.status == "failed":
                    raise RuntimeError(f"transient integration failed at {before:.6g} s")
                steps += 1
                dense = solver.dense_output()
                here = run.at(solver.y[: run.count])
                if here.clearance <= options.contact_clearance:
                    moment = _contact_time(run, dense, before, solver.t, options.contact_clearance)
                    while pending < len(stamps) and stamps[pending] < moment:
                        frames.append((float(stamps[pending]), dense(stamps[pending])))
                        current.advance()
                        pending += 1
                    final = dense(moment)
                    frames.append((moment, final))
                    contact = run.at(final[: run.count]).closest or (-1, -1)
                    break
                while pending < len(stamps) and stamps[pending] <= solver.t:
                    frames.append((float(stamps[pending]), dense(stamps[pending])))
                    current.advance()
                    pending += 1
            if contact is not None:
                break
            y = solver.y
        current.count(steps=steps, evaluations=run.evaluations)
    return _assemble(run, system.dofs, frames, contact, steps, stored, structure)


def _solver(
    run: _Run, time: float, y: Array, end: float, tolerance: float, scale: Array, limit: float
) -> scipy.integrate.Radau:
    return scipy.integrate.Radau(
        run.fun,
        time,
        y,
        end,
        rtol=tolerance,
        atol=scale,  # type: ignore[arg-type]
        jac=run.jac,
        max_step=limit,
        first_step=None if np.isinf(limit) else limit,
    )


def _contact_time(
    run: _Run, dense: Callable[[float], Array], before: float, after: float, threshold: float
) -> float:
    def margin(time: float) -> float:
        try:
            return run.at(dense(time)[: run.count]).clearance - threshold
        except _ClosedError:
            return -threshold

    if margin(before) <= 0:
        return before
    found: float = scipy.optimize.brentq(
        margin, before, after, xtol=1e-12 * (after - before + 1e-300)
    )
    return found


def _rest_charge(run: _Run, place: Array, time: float) -> Array:
    """Charge on each node behind resistance with the sources held without current,
    equation (147)."""
    here = run.at(place)
    joined = run.tie.T @ here.matrix @ run.tie
    nodes = np.zeros(len(joined))
    emf = run.emf(time)
    nodes[run.ideal] = emf[run.held]
    driven = len(run.conductance)
    dynamic = run.dynamic
    if dynamic.size:
        source_emf = np.zeros(driven)
        source_emf[: len(run.sources)] = emf
        right = run.feed[dynamic] * source_emf[dynamic] - run.total[np.ix_(dynamic, run.ideal)] @ (
            nodes[run.ideal]
        )
        nodes[dynamic] = scipy.linalg.solve(run.total[np.ix_(dynamic, dynamic)], right)
    hidden = np.arange(driven, len(joined))
    if hidden.size:
        nodes[hidden] = scipy.linalg.solve(
            joined[np.ix_(hidden, hidden)],
            run.floating - joined[np.ix_(hidden, np.arange(driven))] @ nodes[:driven],
        )
    charge: Array = (joined @ nodes)[dynamic]
    return charge


def _scales(
    structure: Structure,
    run: _Run,
    dofs: tuple[Dof, ...],
    natural: Array,
    stamps: Array,
    tolerance: float,
) -> Array:
    """Absolute tolerance of each state, equation (146)."""
    gap = structure.stack.box.thickness
    lengths = []
    for body, name in (dof for dof, keep in zip(dofs, run.kept, strict=True) if keep):
        if name.startswith("u"):
            lengths.append(gap)
        else:
            area = structure.bodies[body].device.area or structure.bodies[body].handle.area
            lengths.append(gap / max(np.sqrt(area), gap))
    place = np.array(lengths)
    probe = [*stamps, *(p for s in run.sources for p in getattr(s.waveform, "breakpoints", ()))]
    volts = max(1.0, max((float(np.max(abs(run.emf(t)))) for t in probe), default=0.0))
    joined = run.tie.T @ run.model.matrix @ run.tie
    held = np.abs(np.diag(joined))[run.dynamic] * volts
    energy = 0.5 * volts**2 * float(np.sum(np.abs(np.diag(joined)))) or 1e-30
    scale: Array = tolerance * np.concatenate(
        [place, place * float(natural[0]), held, [energy] * 3]
    )
    return scale


def _origin(model: CapacitanceModel, start: Equilibrium | ArrayLike | None, kept: Array) -> Array:
    if start is None:
        return np.zeros(len(model.dofs))
    if isinstance(start, Equilibrium):
        if start.model.rest.structure is not model.rest.structure:
            raise ValueError("start is an equilibrium of another structure")
        return start.displacement
    place = np.asarray(start, dtype=float)
    if place.shape != (len(model.dofs),):
        raise ValueError(
            f"start must hold one value for each of {len(model.dofs)} degrees of freedom"
        )
    return np.where(kept, place, 0.0)


def _assemble(
    run: _Run,
    dofs: tuple[Dof, ...],
    frames: list[tuple[float, Array]],
    contact: tuple[int, int] | None,
    steps: int,
    stored: Mapping[int, float],
    structure: Structure,
) -> Transient:
    n = run.count
    columns: dict[str, list[Any]] = {key: [] for key in ("place", "speed", "volts", "charge")}
    extra: dict[str, list[float]] = {
        key: [] for key in ("kinetic", "strain", "field", "damped", "lost", "supplied", "gap")
    }
    currents = []
    for time, y in frames:
        now = run.state(time, y)
        place = np.zeros(len(dofs))
        speed = np.zeros(len(dofs))
        place[run.kept], speed[run.kept] = y[:n], y[n : 2 * n]
        columns["place"].append(place)
        columns["speed"].append(speed)
        columns["volts"].append(now.electrodes)
        columns["charge"].append(now.here.matrix @ now.electrodes)
        currents.append(now.current)
        extra["kinetic"].append(0.5 * float(y[n : 2 * n] @ run.mass @ y[n : 2 * n]))
        extra["strain"].append(0.5 * float(y[:n] @ run.spring @ y[:n]))
        extra["field"].append(0.5 * float(now.electrodes @ now.here.matrix @ now.electrodes))
        extra["damped"].append(float(y[-3]))
        extra["lost"].append(float(y[-2]))
        extra["supplied"].append(float(y[-1]))
        extra["gap"].append(now.here.clearance)
    last = run.state(*frames[-1]).here
    return Transient(
        times=np.array([time for time, _ in frames]),
        displacement=np.array(columns["place"]),
        velocity=np.array(columns["speed"]),
        potentials=np.array(columns["volts"]),
        charge=np.array(columns["charge"]),
        current=np.array(currents).reshape(len(frames), len(run.sources)),
        clearance=np.array(extra["gap"]),
        energy=Energy(
            kinetic=np.array(extra["kinetic"]),
            strain=np.array(extra["strain"]),
            electrostatic=np.array(extra["field"]),
            damped=np.array(extra["damped"]),
            resistive=np.array(extra["lost"]),
            supplied=np.array(extra["supplied"]),
        ),
        contact=contact,
        dofs=dofs,
        loose=tuple(
            sorted({body for (body, _), keep in zip(dofs, run.kept, strict=True) if not keep})
        ),
        floating=tuple(sorted(stored)),
        steps=steps,
        evaluations=run.evaluations,
        model=last,
    )


def _value(waveform: float | Callable[[float], float], time: float) -> float:
    if callable(waveform):
        return float(waveform(time))
    return float(waveform)


def _slope(waveform: float | Callable[[float], float], time: float, duration: float) -> float:
    if not callable(waveform):
        return 0.0
    derivative = getattr(waveform, "derivative", None)
    if derivative is not None:
        return float(derivative(time))
    step = 1e-7 * max(duration, abs(time), 1e-300)
    return (float(waveform(time + step)) - float(waveform(time - step))) / (2 * step)

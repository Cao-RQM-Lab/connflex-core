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
from connflex_core.electrical.drivers import Capacitor, CurrentSource, Resistor, Source, Switch
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
        electrostatic: Energy in the field between electrodes and in capacitors of the
            circuit.
        damped: Energy that damping has taken since the start.
        resistive: Energy that resistance of links, sources, and resistors of the circuit
            has taken since the start.
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
            the order of the sources: for a current source, its waveform.
        clearance: Clearance of equation (126) at each time.
        energy: Energy ledger at each time.
        contact: Electrodes across the gap whose closing to ``contact_clearance`` ended
            the run. None when the run reached its last time. ``(-1, -1)`` when the
            integrator could not step on from a state with no gap narrower than at
            rest, as from a start beyond pull-in: the last frame is that state, with
            its clearance still above ``contact_clearance``.
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
    currents: Sequence[CurrentSource]
    inject: Array
    order: list[tuple[bool, int]]
    tie: Array
    ideal: NDArray[np.intp]
    dynamic: NDArray[np.intp]
    unknown: NDArray[np.intp]
    floating: Array
    links: Array
    resistors: list[tuple[Array, Switch | None]]
    feeds: Array
    lumped: Array
    kept: NDArray[np.bool_]
    mass: Array
    spring: Array
    damping: Array
    applied: Array
    duration: float
    evaluations: int = 0
    cache: dict[bytes, CapacitanceModel] = field(default_factory=dict)
    workers: int | None = None

    def __post_init__(self) -> None:
        self.factor = scipy.linalg.cho_factor(self.mass)
        self.count = int(self.kept.sum())
        self.dofs = len(self.kept)
        self.held = np.array([source.resistance == 0 for source in self.sources], dtype=bool)
        self.slowed = scipy.linalg.cho_solve(self.factor, self.damping)
        self.switch_to(-np.inf)

    def switch_to(self, time: float) -> None:
        """Set the conductance of the network with switches as they are at ``time``."""
        self.conductance = self.links + sum(
            (stamp for stamp, switch in self.resistors if switch is None or switch.closed_at(time)),
            np.zeros_like(self.links),
        )
        self.feed = self.feeds.copy()
        for k, source in enumerate(self.sources):
            if source.switch is not None and not source.switch.closed_at(time):
                self.feed[k] = 0.0
        self.total = self.conductance + np.diag(self.feed)

    def joined(self, here: CapacitanceModel) -> Array:
        """Capacitance between nodes: of the structure and of the circuit, equation (142)."""
        joined: Array = self.tie.T @ here.matrix @ self.tie + self.lumped
        return joined

    def pushed(self, time: float) -> Array:
        """Current of each current source at ``time``."""
        return np.array([_value(source.waveform, time) for source in self.currents], dtype=float)

    def injected(self, pushed: Array) -> Array:
        """Current that current sources of ``pushed`` drive into each node."""
        injected: Array = self.inject @ pushed if len(pushed) else np.zeros(len(self.links))
        return injected

    def at(self, place: Array) -> CapacitanceModel:
        key = place.tobytes()
        if key not in self.cache:
            full = np.zeros(self.dofs)
            full[self.kept] = place
            try:
                self.cache[key] = self.model.at(full, workers=self.workers)
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
        joined = self.joined(here)
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
        forced = self.pushed(time)
        injected = self.injected(forced)
        flowing = -self.total @ nodes[:driven] + self.feed * source_emf + injected
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
            drift[ideal]
            + joined[ideal] @ rising
            + (self.conductance @ nodes[:driven])[ideal]
            - injected[ideal]
        )
        driving = np.zeros(len(self.sources))
        driving[self.held] = held_current
        resistive = ~self.held
        driving[resistive] = (self.feed * (source_emf - nodes[:driven]))[: len(self.sources)][
            resistive
        ]
        current = np.array(
            [driving[k] if voltage else forced[k] for voltage, k in self.order], dtype=float
        )

        damped = float(speed @ self.damping @ speed)
        across = nodes[:driven]
        lost = float(across @ self.conductance @ across) + float(
            self.feed @ (source_emf - across) ** 2
        )
        supplied = (
            float(emf @ driving) + float(injected @ nodes[:driven]) + float(self.applied @ speed)
        )
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
    sources: Sequence[Source | CurrentSource],
    times: ArrayLike,
    options: TransientOptions,
    *,
    circuit: Sequence[Resistor | Capacitor] = (),
    charges: Mapping[int, float] | None = None,
    loads: Mapping[int, Sequence[float]] | None = None,
    start: Equilibrium | ArrayLike | None = None,
    workers: int | None = None,
    progress: Progress | None = None,
) -> Transient:
    """Motion, potentials, and charges of ``structure`` over ``times`` under ``sources``.

    ``model`` is the capacitance model of ``structure`` at rest. Each voltage or current
    source drives one electrode; electrodes of its net that links with resistance set
    apart from it are reached through them. ``circuit`` holds resistors and capacitors
    between electrodes or to ground. ``charges`` maps a floating net to the charge it
    keeps. A net that no source, element of the circuit, or charge reaches is held at
    zero volts on every electrode. ``loads`` are constant forces on bodies, as for
    ``static``. The run starts at ``times[0]`` from rest or from the displacement of
    ``start``, at zero velocity, with the potentials that the sources hold without
    current at that time, equation (147), and ends at ``times[-1]`` or at contact. A
    start at which bodies touch is refused with the model's error.
    Integration restarts where a waveform or switch changes. ``workers`` threads
    evaluate the model, as many as at rest by default. Stages assemble and transient
    log through the ``connflex_core`` logger and report to ``progress``, which advances
    once per entry of ``times`` of the result, so a contact event leaves the stage
    short of its total.
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
        if not isinstance(source, Source | CurrentSource):
            raise TypeError(
                f"a source is a Source or a CurrentSource, not {type(source).__name__} {source!r}"
            )
    parts: Sequence[object] = circuit
    for part in parts:
        if not isinstance(part, Resistor | Capacitor):
            raise TypeError(
                f"an element of the circuit is a Resistor or Capacitor, not "
                f"{type(part).__name__} {part!r}"
            )
    voltage = [source for source in sources if isinstance(source, Source)]
    currents = [source for source in sources if isinstance(source, CurrentSource)]
    order, counts = [], [0, 0]
    for source in sources:
        kind = isinstance(source, Source)
        order.append((kind, counts[kind]))
        counts[kind] += 1
    places = electrode_indices(
        structure, [source.electrode for source in voltage], stored, "source"
    )
    feeding = electrode_indices(
        structure, [source.electrode for source in currents], stored, "current source"
    )
    ends = [
        electrode_indices(
            structure,
            [part.electrode] if part.other is None else [part.electrode, part.other],
            stored,
            "end of a circuit element",
        )
        for part in circuit
    ]
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

    reached = [*feeding, *(end for pair in ends for end in pair)]
    nodes, hidden, conductance = network_nodes(structure, stored, places, "source", reached)
    tie = np.zeros((len(structure.electrodes), len(nodes) + len(hidden)))
    for column, node in enumerate([*nodes, *hidden]):
        tie[list(node), column] = 1.0
    index = {electrode: k for k, node in enumerate(nodes) for electrode in node}
    resistors, lumped = [], np.zeros((len(tie[0]), len(tie[0])))
    for part, pair in zip(circuit, ends, strict=True):
        at = [index[end] for end in pair]
        if len(at) == 2 and at[0] == at[1]:
            raise ValueError(
                f"both ends of a circuit element, electrodes {pair[0]} and {pair[1]}, are one "
                f"node: links without resistance join them"
            )
        if isinstance(part, Resistor):
            resistors.append((_stamp(len(nodes), at, 1 / part.resistance), part.switch))
        else:
            lumped += _stamp(len(lumped), at, part.capacitance)
    inject = np.zeros((len(nodes), len(currents)))
    for k, electrode in enumerate(feeding):
        inject[index[electrode], k] = 1.0
    held = np.array([source.resistance == 0 for source in voltage], dtype=bool)
    ideal = np.flatnonzero(held).astype(np.intp)
    dynamic = np.setdiff1d(np.arange(len(nodes)), ideal).astype(np.intp)
    unknown = np.concatenate([dynamic, np.arange(len(nodes), len(nodes) + len(hidden))])
    feed = np.zeros(len(nodes))
    for k, source in enumerate(voltage):
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
        workers=workers,
        sources=voltage,
        currents=currents,
        inject=inject,
        order=order,
        tie=tie,
        ideal=ideal,
        dynamic=dynamic,
        unknown=unknown.astype(np.intp),
        floating=floating,
        links=conductance,
        resistors=resistors,
        feeds=feed,
        lumped=lumped,
        kept=kept,
        mass=mass,
        spring=spring,
        damping=damping,
        applied=applied[kept],
        duration=float(stamps[-1] - stamps[0]),
    )
    start_place = origin[kept]
    try:
        run.at(start_place)
    except _ClosedError as exc:  # the model says what closes at this start
        raise ValueError(str(exc.__cause__)) from exc.__cause__
    y = np.concatenate(
        [start_place, np.zeros(run.count), _rest_charge(run, start_place, stamps[0]), np.zeros(3)]
    )
    scale = _scales(structure, run, system.dofs, natural, stamps, options.tolerance)
    switches = [part.switch for part in [*voltage, *circuit] if isinstance(part, Source | Resistor)]
    changes = [
        *(getattr(source.waveform, "breakpoints", ()) for source in sources),
        *(switch.breakpoints for switch in switches if switch is not None),
    ]
    breaks = sorted(
        {float(point) for points in changes for point in points if stamps[0] < point < stamps[-1]}
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
            run.switch_to(left)
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
                        current.advance()
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
                    current.advance()
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
    equation (147).

    Nodes that resistance ties to a source or ground take the potentials that carry no
    current. Nodes that resistance joins to neither share one potential per group and
    hold no charge together, and floating nets keep their charge.
    """
    run.switch_to(time)
    here = run.at(place)
    joined = run.joined(here)
    nodes = np.zeros(len(joined))
    emf = run.emf(time)
    nodes[run.ideal] = emf[run.held]
    driven, dynamic = len(run.conductance), run.dynamic
    source_emf = np.zeros(driven)
    source_emf[: len(run.sources)] = emf
    within = run.total[np.ix_(dynamic, dynamic)]
    leaks = within.sum(axis=1) > 1e-9 * float(np.max(np.abs(np.diag(run.total)), initial=0.0))
    graph: nx.Graph[int] = nx.Graph()
    graph.add_nodes_from(range(len(dynamic)))
    graph.add_edges_from((int(i), int(j)) for i, j in np.argwhere(np.triu(within != 0, 1)))
    groups = [sorted(group) for group in nx.connected_components(graph)]
    tied = np.array(sorted(k for g in groups if leaks[g].any() for k in g), dtype=np.intp)
    if tied.size:
        tied_nodes = dynamic[tied]
        right = (
            run.feed[tied_nodes] * source_emf[tied_nodes]
            + run.injected(run.pushed(time))[tied_nodes]
            - run.total[np.ix_(tied_nodes, run.ideal)] @ nodes[run.ideal]
        )
        nodes[tied_nodes] = scipy.linalg.solve(run.total[np.ix_(tied_nodes, tied_nodes)], right)
    loose = [dynamic[g] for g in groups if not leaks[g].any()]
    hidden = np.arange(driven, len(joined))
    if loose or hidden.size:
        free = np.concatenate([*loose, hidden]).astype(np.intp)
        share = np.zeros((len(free), len(loose) + len(hidden)))
        row = 0
        for column, group in enumerate(loose):
            share[row : row + len(group), column] = 1.0
            row += len(group)
        share[row:, len(loose) :] = np.eye(len(hidden))
        known = np.setdiff1d(np.arange(len(joined)), free)
        target = np.concatenate([np.zeros(row), run.floating])
        right = share.T @ (target - joined[np.ix_(free, known)] @ nodes[known])
        nodes[free] = share @ scipy.linalg.solve(
            share.T @ joined[np.ix_(free, free)] @ share, right
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
    volts = max(1.0, *(float(np.max(abs(run.emf(t)), initial=0.0)) for t in probe))
    joined = run.joined(run.model)
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
        run.switch_to(time)
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
        extra["field"].append(0.5 * float(now.potentials @ now.joined @ now.potentials))
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
        current=np.array(currents).reshape(len(frames), len(run.order)),
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


def _stamp(size: int, ends: list[int], value: float) -> Array:
    """Matrix of a branch of ``value`` between two nodes, or from one node to ground."""
    stamp = np.zeros((size, size))
    for k in ends:
        stamp[k, k] += value
    if len(ends) == 2:
        stamp[ends[0], ends[1]] -= value
        stamp[ends[1], ends[0]] -= value
    return stamp


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

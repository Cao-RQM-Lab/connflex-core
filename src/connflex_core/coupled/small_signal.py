"""Small-signal response about a static equilibrium."""

from __future__ import annotations

import warnings
from collections.abc import Sequence
from dataclasses import dataclass, field
from itertools import pairwise
from typing import TYPE_CHECKING

import networkx as nx
import numpy as np
import scipy.linalg
from numpy.typing import ArrayLike, NDArray

from connflex_core.coupled.static import Equilibrium
from connflex_core.electrical.electrodes import Electrode
from connflex_core.mechanics.assembly import Dof, assemble
from connflex_core.params import SmallSignalOptions
from connflex_core.progress import Progress, stage

if TYPE_CHECKING:
    from connflex_core.model import Structure

Array = NDArray[np.float64]
Complex = NDArray[np.complex128]


class LinearizationWarning(UserWarning):
    """A drive reaches beyond the range where the small-signal model holds."""


@dataclass(frozen=True)
class Linearity:
    """How far a drive at one frequency strains the small-signal model, equation (137).

    Attributes:
        voltage: Largest ratio, over pairs of electrodes whose capacitance changes with
            motion, of the swing of the potential between them to their potential
            difference at the equilibrium. Infinite when a pair without bias swings.
        gap: Fraction by which the narrowest gap between two bodies closes at the peak
            of the motion, against its width at the equilibrium. 1 when the motion
            would close it.
        peak: Displacement of every degree of freedom at the moment of the cycle when
            the motion stores most strain energy in the beams.
    """

    voltage: float
    gap: float
    peak: Array


@dataclass(frozen=True, eq=False)
class _Network:
    """Nodes of the network and the matrices that reconstruct floating potentials."""

    tie: Array
    floating: Array
    within: Array
    across: Array
    pull: Array
    kept: NDArray[np.bool_]
    spring: Array


@dataclass(frozen=True, eq=False)
class SmallSignal:
    """Linear response about an equilibrium to small potentials at ports, equations (129)
    to (137).

    Attributes:
        ports: Electrode of each port, in the order given. A source holds each at the
            potential of its net in the equilibrium plus the signal.
        nodes: Indices of the electrodes of each node of the network: the ports first,
            then the electrodes of their nets that links with resistance set apart from
            them. Links without resistance join electrodes into one node.
        frequencies_hz: Resonant frequencies at the equilibrium, ascending, equation (132).
        modes: One column per mode over ``dofs``, mass-normalized, zero on loose bodies.
        dofs: (body index, degree of freedom) of each row of ``modes``.
        transduction: Mode by node: generalized force on each mode per volt on each
            node, equation (132). Zero below the rounding of the mode shapes.
        capacitance: Node by node capacitance, floating nets eliminated, equation (131).
        conductance: Node by node conductance of links with resistance, equation (130).
        options: Options it was built with, the quality factor among them.
        equilibrium: Equilibrium the model is linearized about.
    """

    ports: tuple[Electrode, ...]
    nodes: tuple[tuple[int, ...], ...]
    frequencies_hz: Array
    modes: Array
    dofs: tuple[Dof, ...]
    transduction: Array
    capacitance: Array
    conductance: Array
    options: SmallSignalOptions
    equilibrium: Equilibrium
    _network: _Network = field(repr=False)

    def admittance(self, frequencies: ArrayLike) -> Complex:
        """Port by port admittance at each frequency in hertz, equations (134) and (135).

        Entry (f, p, r) is the current a source drives into port p per volt at port r,
        with every other port held at its bias.
        """
        omega = _angular(frequencies)
        nodes = self._nodes(omega)
        count = len(self.ports)
        inner = slice(count, None)
        outer = slice(None, count)
        if nodes.shape[1] == count:
            return nodes
        solved = np.linalg.solve(nodes[:, inner, inner], nodes[:, inner, outer])
        reduced: Complex = nodes[:, outer, outer] - nodes[:, outer, inner] @ solved
        return reduced

    def displacement(self, frequencies: ArrayLike, drive: ArrayLike) -> Complex:
        """Complex amplitude of every degree of freedom at each frequency in hertz for
        ``drive``, the complex amplitude in volts at each port, equation (133)."""
        omega = _angular(frequencies)
        swing = self._swing(omega, drive)
        response = self._response(omega)
        motion: Complex = np.einsum(
            "dn,fn,ns,fs->fd", self.modes, response, self.transduction, swing
        )
        return motion

    def motional(self) -> tuple[Array, Array, Array]:
        """Inductance, capacitance, and resistance of the series branch each mode adds
        between each node and ground with every other node held, equation (136); mode
        by node. A node that a mode does not couple to has infinite inductance."""
        coupling = self.transduction**2
        natural = (2 * np.pi * self.frequencies_hz)[:, None]
        with np.errstate(divide="ignore"):
            inductance = 1 / coupling
        capacitance = coupling / natural**2
        resistance = natural * inductance / self.options.quality_factor
        return inductance, capacitance, resistance

    def linearity(self, drive: ArrayLike, frequency: float) -> Linearity:
        """Voltage and gap ratios of equation (137) for ``drive`` at ``frequency`` in hertz.

        Warns with ``LinearizationWarning`` when either exceeds ``linear_tolerance``. The
        gap ratio evaluates the capacitance model twice, at either extreme of the motion.
        """
        omega = _angular([frequency])
        swing = self._swing(omega, drive)[0]
        response = self._response(omega)[0]
        motion = np.einsum("dn,n,ns,s->d", self.modes, response, self.transduction, swing)
        network = self._network
        signal = network.tie[:, : len(self.nodes)]
        potentials = signal @ swing
        if network.floating.size:
            held = network.across @ swing + network.pull @ motion
            potentials += network.floating @ scipy.linalg.solve(network.within, -held)
        voltage = _voltage_ratio(self.equilibrium, potentials, network.kept)
        peak = _peak(motion, network.spring)
        model = self.equilibrium.model
        gap = 0.0
        for sign in (1.0, -1.0):
            try:
                moved = model.at(model.displacement + sign * peak).clearance
            except ValueError:
                gap = 1.0
                break
            gap = max(gap, 1.0 - moved / model.clearance)
        tolerance = self.options.linear_tolerance
        for name, value in (("voltage", voltage), ("gap", gap)):
            if value > tolerance:
                warnings.warn(
                    f"{name} ratio {value:.3g} exceeds linear_tolerance {tolerance:g}: "
                    f"the drive reaches beyond the small-signal range",
                    LinearizationWarning,
                    stacklevel=2,
                )
        return Linearity(voltage=voltage, gap=gap, peak=peak)

    def _response(self, omega: Array) -> Complex:
        natural = 2 * np.pi * self.frequencies_hz
        denominator = (
            natural**2
            - omega[:, None] ** 2
            + 1j * omega[:, None] * natural / self.options.quality_factor
        )
        response: Complex = np.reciprocal(denominator.astype(complex))
        return response

    def _nodes(self, omega: Array) -> Complex:
        """Node by node admittance at each angular frequency, equation (134)."""
        coupled = np.einsum(
            "ns,fn,nt->fst", self.transduction, self._response(omega), self.transduction
        )
        nodes: Complex = self.conductance + 1j * omega[:, None, None] * (self.capacitance + coupled)
        return nodes

    def _swing(self, omega: Array, drive: ArrayLike) -> Complex:
        """Potential amplitude of every node for ``drive`` at the ports, equation (135)."""
        given = np.asarray(drive, dtype=complex)
        count = len(self.ports)
        if given.shape != (count,):
            raise ValueError(f"drive must hold one amplitude for each of {count} ports")
        nodes = self._nodes(omega)
        swing = np.zeros((len(omega), len(self.nodes)), dtype=complex)
        swing[:, :count] = given
        if len(self.nodes) > count:
            pushed = nodes[:, count:, :count] @ given
            swing[:, count:] = -np.linalg.solve(nodes[:, count:, count:], pushed[..., None])[..., 0]
        return swing


def small_signal(
    structure: Structure,
    equilibrium: Equilibrium,
    ports: Sequence[Electrode],
    options: SmallSignalOptions,
    *,
    progress: Progress | None = None,
) -> SmallSignal:
    """Small-signal model of ``structure`` about ``equilibrium`` with sources at ``ports``.

    ``ports`` are electrodes of ``structure``, as ``Structure.electrode_at`` returns
    them, each held by a source at the potential its net has in the equilibrium plus the
    small signal; an integer is refused, so a net cannot pass for one. Electrodes of a
    port's net that links with resistance set apart from it are inner nodes. A net with
    a potential and no port is held at every electrode. A floating net keeps its charge,
    is one node, and cannot carry a port. Stages assemble and small signal log through
    the ``connflex_core`` logger and report to ``progress``.
    """
    model = equilibrium.model
    if model.rest.structure is not structure:
        raise ValueError("equilibrium is of another structure")
    if not equilibrium.stable:
        raise ValueError("small signal needs a stable equilibrium; this one does not hold")
    count = len(structure.electrodes)
    if not ports:
        raise ValueError("small signal needs at least one port")
    given: Sequence[object] = ports
    places = []
    for port in given:
        if not isinstance(port, Electrode):
            raise TypeError(
                f"a port is an Electrode, such as structure.electrode_at(point, layer), not "
                f"{type(port).__name__} {port!r}: electrodes and nets are numbered apart"
            )
        if not (0 <= port.index < count and structure.electrodes[port.index] is port):
            raise ValueError(f"electrode {port.index} is of another structure")
        if port.net in equilibrium.floating:
            raise ValueError(
                f"electrode {port.index} is on floating net {port.net}, which no source holds"
            )
        places.append(port.index)
    with stage("assemble", progress) as current:
        system = assemble(structure)
        current.count(dofs=len(system.dofs), beams=len(structure.beams))
    with stage("small signal", progress) as current:
        nodes, floating, conductance = _network_nodes(structure, equilibrium, places)
        tie = np.zeros((count, len(nodes) + len(floating)))
        for column, node in enumerate([*nodes, *floating]):
            tie[list(node), column] = 1.0
        pull = tie.T @ np.einsum("abk,b->ak", model.gradient, equilibrium.potentials)
        joined = tie.T @ model.matrix @ tie
        signal = slice(len(nodes))
        hidden = slice(len(nodes), None)
        if floating:
            factor = scipy.linalg.cho_factor(joined[hidden, hidden])
            capacitance = joined[signal, signal] - joined[signal, hidden] @ scipy.linalg.cho_solve(
                factor, joined[hidden, signal]
            )
            coupling = pull[signal] - joined[signal, hidden] @ scipy.linalg.cho_solve(
                factor, pull[hidden]
            )
        else:
            capacitance, coupling = joined, pull
        kept = np.array([body not in equilibrium.loose for body, _ in system.dofs])
        both = np.ix_(kept, kept)
        squared, shapes = scipy.linalg.eigh(equilibrium.stiffness[both], system.mass[both])
        modes = np.zeros((len(system.dofs), len(squared)))
        modes[kept] = shapes
        transduction = (coupling @ modes).T
        floor = kept.sum() * np.finfo(float).eps * np.linalg.norm(transduction, axis=0)
        transduction[abs(transduction) <= floor] = 0.0
        current.count(modes=len(squared), nodes=len(nodes))
    return SmallSignal(
        ports=tuple(ports),
        nodes=tuple(nodes),
        frequencies_hz=np.sqrt(squared) / (2 * np.pi),
        modes=modes,
        dofs=system.dofs,
        transduction=transduction,
        capacitance=(capacitance + capacitance.T) / 2,
        conductance=conductance,
        options=options,
        equilibrium=equilibrium,
        _network=_Network(
            tie=tie,
            floating=tie[:, hidden],
            within=joined[hidden, hidden],
            across=joined[hidden, signal],
            pull=pull[hidden],
            kept=kept,
            spring=system.stiffness,
        ),
    )


def _network_nodes(
    structure: Structure, equilibrium: Equilibrium, ports: Sequence[int]
) -> tuple[list[tuple[int, ...]], list[tuple[int, ...]], Array]:
    """Signal nodes with ports first, floating nodes, and conductance between signal
    nodes, equation (130)."""
    net = equilibrium.model.net
    floating = set(equilibrium.floating)
    sourced = {int(net[port]) for port in ports}
    shorts: nx.Graph[int] = nx.Graph()
    shorts.add_nodes_from(range(len(structure.electrodes)))
    for link in structure.links:
        if not link.resistance:
            shorts.add_edge(*link.electrodes)
    for each in floating:
        members = np.flatnonzero(net == each).tolist()
        shorts.add_edges_from(pairwise(members))
    groups = {a: tuple(sorted(group)) for group in nx.connected_components(shorts) for a in group}
    nodes: list[tuple[int, ...]] = []
    for port in ports:
        if groups[port] in nodes:
            raise ValueError(
                f"ports at electrodes {nodes[nodes.index(groups[port])][0]} and {port} are one "
                f"node: links without resistance join them"
            )
        nodes.append(groups[port])
    rest = sorted(
        {groups[a] for a in range(len(net)) if net[a] in sourced and groups[a] not in nodes}
    )
    nodes.extend(rest)
    hidden = sorted({groups[a] for a in range(len(net)) if net[a] in floating})
    index = {a: k for k, node in enumerate(nodes) for a in node}
    conductance = np.zeros((len(nodes), len(nodes)))
    for link in structure.links:
        a, b = link.electrodes
        if link.resistance and a in index and index[a] != index[b]:
            pair = [index[a], index[b]]
            conductance[np.ix_(pair, pair)] += np.array([[1, -1], [-1, 1]]) / link.resistance
    return nodes, hidden, conductance


def _angular(frequencies: ArrayLike) -> Array:
    hertz = np.atleast_1d(np.asarray(frequencies, dtype=float))
    if hertz.ndim != 1 or np.any(hertz < 0) or not np.all(np.isfinite(hertz)):
        raise ValueError("frequencies must be a sequence of finite values of zero or more hertz")
    angular: Array = 2 * np.pi * hertz
    return angular


def _voltage_ratio(equilibrium: Equilibrium, swing: Complex, kept: NDArray[np.bool_]) -> float:
    """Largest swing across a pair whose capacitance moves, against its bias."""
    gradient = equilibrium.model.gradient[:, :, kept]
    pairs = np.argwhere(np.triu(np.any(gradient != 0, axis=2), 1))
    if not len(pairs):
        return 0.0
    a, b = pairs.T
    across = abs(swing[a] - swing[b])
    bias = abs(equilibrium.potentials[a] - equilibrium.potentials[b])
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(across > 0, across / bias, 0.0)
    return float(np.max(ratio))


def _peak(motion: Complex, spring: Array) -> Array:
    """Real displacement over a cycle of ``motion`` that stores most strain energy."""
    parts = np.stack([motion.real, -motion.imag], axis=1)
    gram = parts.T @ spring @ parts
    _, vectors = np.linalg.eigh(gram)
    peak: Array = parts @ vectors[:, -1]
    return peak

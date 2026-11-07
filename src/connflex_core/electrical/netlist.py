"""SPICE netlist of a small-signal model."""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from connflex_core.coupled.small_signal import SmallSignal


def spice_netlist(signal: SmallSignal, name: str, modes: Sequence[int] | None) -> str:
    """Subcircuit ``name`` with pins ``port0``, ``port1``, ... in the order of the ports.

    Capacitance and conductance between nodes become capacitors and resistors, and
    each mode a series branch that sources couple to the nodes, equation (138).
    Ground is node 0, the potential of nets held without a port.
    """
    total = len(signal.frequencies_hz)
    chosen = list(range(total)) if modes is None else [int(mode) for mode in modes]
    for mode in chosen:
        if not 0 <= mode < total:
            raise ValueError(f"mode {mode} does not exist; the model has {total}")
    written = [mode for mode in chosen if np.max(abs(signal.transduction[mode]))]
    ports = len(signal.ports)
    pins = [f"port{k}" for k in range(ports)]
    names = pins + [f"n{k}" for k in range(ports, len(signal.nodes))]
    quality = signal.options.quality_factor
    lines = [
        f"* small-signal model about an equilibrium: {ports} ports, {len(names)} nodes, "
        f"{len(written)} of {total} modes, quality factor {quality:g}",
        *(f"* {pin}: electrode {port.index}" for pin, port in zip(pins, signal.ports, strict=True)),
        f".subckt {name} {' '.join(pins)}",
    ]
    capacitance, conductance = signal.capacitance, signal.conductance
    for i in range(len(names)):
        if grounded := capacitance[i].sum():
            lines.append(f"Cg{i} {names[i]} 0 {_number(grounded)}")
        for j in range(i + 1, len(names)):
            if capacitance[i, j]:
                lines.append(f"Cc{i}_{j} {names[i]} {names[j]} {_number(-capacitance[i, j])}")
            if conductance[i, j]:
                lines.append(f"Rn{i}_{j} {names[i]} {names[j]} {_number(-1 / conductance[i, j])}")
    for mode in written:
        coupling = signal.transduction[mode]
        scale = float(np.max(abs(coupling)))
        natural = 2 * np.pi * float(signal.frequencies_hz[mode])
        lines.append(f"* mode {mode} at {signal.frequencies_hz[mode]:.9g} Hz")
        previous = "0"
        coupled = [node for node in range(len(names)) if coupling[node]]
        for node in coupled:
            gain = _number(coupling[node] / scale)
            lines.append(f"Em{mode}_{node} m{mode}e{node} {previous} {names[node]} 0 {gain}")
            previous = f"m{mode}e{node}"
        lines += [
            f"Lm{mode} {previous} m{mode}l {_number(1 / scale**2)}",
            f"Rm{mode} m{mode}l m{mode}r {_number(natural / (quality * scale**2))}",
            f"Cm{mode} m{mode}r m{mode}c {_number(scale**2 / natural**2)}",
            f"Vm{mode} m{mode}c 0 0",
        ]
        for node in coupled:
            gain = _number(coupling[node] / scale)
            lines.append(f"Fm{mode}_{node} {names[node]} 0 Vm{mode} {gain}")
    lines.append(f".ends {name}")
    return "\n".join(lines) + "\n"


def _number(value: float) -> str:
    return f"{float(value):.12g}"

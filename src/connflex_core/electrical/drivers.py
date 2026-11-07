"""Sources and circuit elements that drive electrodes over time, and their waveforms."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np

from connflex_core._validation import require_ascending, require_nonnegative, require_positive
from connflex_core.electrical.electrodes import Electrode

Waveform = float | Callable[[float], float]


@dataclass(frozen=True)
class PiecewiseLinear:
    """Waveform through points, constant before the first and after the last.

    Attributes:
        times: Times of the points in seconds, strictly ascending.
        values: Value at each time, in volts for a source.
    """

    times: tuple[float, ...]
    values: tuple[float, ...]

    def __init__(self, times: Sequence[float], values: Sequence[float]) -> None:
        stamps = require_ascending("times", times, 1)
        levels = tuple(float(value) for value in values)
        if len(stamps) != len(levels):
            raise ValueError("times and values must be equally long")
        if not np.all(np.isfinite(levels)):
            raise ValueError("values must be finite")
        object.__setattr__(self, "times", stamps)
        object.__setattr__(self, "values", levels)

    def __call__(self, time: float) -> float:
        return float(np.interp(time, self.times, self.values))

    def derivative(self, time: float) -> float:
        """Slope at ``time``, that of the segment that starts there at a point."""
        segment = int(np.searchsorted(self.times, time, side="right")) - 1
        if segment < 0 or segment >= len(self.times) - 1:
            return 0.0
        rise = self.values[segment + 1] - self.values[segment]
        return rise / (self.times[segment + 1] - self.times[segment])

    @property
    def breakpoints(self) -> tuple[float, ...]:
        """Times where the slope changes."""
        return self.times


@dataclass(frozen=True)
class Switch:
    """Times at which a branch with resistance opens or closes, equation (141).

    Attributes:
        times: Times in seconds, strictly ascending, at each of which the branch changes
            from closed to open or back. At one of these times it has changed.
        closed: Whether the branch conducts before the first time.
    """

    times: tuple[float, ...]
    closed: bool = True

    def __init__(self, times: Sequence[float], closed: bool = True) -> None:
        object.__setattr__(self, "times", require_ascending("switch times", times, 0))
        object.__setattr__(self, "closed", bool(closed))

    def closed_at(self, time: float) -> bool:
        """Whether the branch conducts at ``time``."""
        changes = int(np.searchsorted(self.times, time, side="right"))
        return self.closed != bool(changes % 2)

    @property
    def breakpoints(self) -> tuple[float, ...]:
        """Times where the branch changes."""
        return self.times


@dataclass(frozen=True)
class Source:
    """Voltage source that drives one electrode in a transient, equation (141).

    Attributes:
        electrode: Electrode it connects to, as ``Structure.electrode_at`` returns it.
        waveform: Volts against seconds: a number for a constant, a ``PiecewiseLinear``,
            or any continuous callable of time. Integration restarts at a waveform's
            ``breakpoints`` where it has them, and a ``derivative`` gives the slope for
            the energy an ideal source supplies; a callable without one is
            differentiated numerically.
        resistance: Series resistance in ohms. Zero holds the electrode at the waveform.
        switch: Times at which the source connects and disconnects; a switched source
            needs resistance.
    """

    electrode: Electrode
    waveform: Waveform
    resistance: float = 0.0
    switch: Switch | None = None

    def __post_init__(self) -> None:
        require_nonnegative("resistance", self.resistance)
        if self.switch is not None and not self.resistance:
            raise ValueError("a switched source needs resistance: an ideal source holds its node")


@dataclass(frozen=True)
class CurrentSource:
    """Current source that drives one electrode from ground in a transient, equation (141).

    Attributes:
        electrode: Electrode it drives current into.
        waveform: Amperes against seconds, as for ``Source``.
    """

    electrode: Electrode
    waveform: Waveform


@dataclass(frozen=True)
class Resistor:
    """Resistor outside the structure, between two electrodes or one and ground, equation
    (141).

    Attributes:
        electrode: Electrode at one end.
        other: Electrode at the other end, or None for ground at zero volts.
        resistance: Resistance in ohms, positive.
        switch: Times at which it connects and disconnects; always connected without one.
    """

    electrode: Electrode
    other: Electrode | None
    resistance: float
    switch: Switch | None = None

    def __post_init__(self) -> None:
        require_positive("resistance", self.resistance)


@dataclass(frozen=True)
class Capacitor:
    """Capacitor outside the structure, between two electrodes or one and ground,
    equation (142).

    Attributes:
        electrode: Electrode at one end.
        other: Electrode at the other end, or None for ground at zero volts.
        capacitance: Capacitance in farads, positive.
    """

    electrode: Electrode
    other: Electrode | None
    capacitance: float

    def __post_init__(self) -> None:
        require_positive("capacitance", self.capacitance)

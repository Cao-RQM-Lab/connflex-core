"""Sources that drive electrodes over time, and their waveforms."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np

from connflex_core._validation import require_nonnegative
from connflex_core.electrical.electrodes import Electrode


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
        stamps = tuple(float(time) for time in times)
        levels = tuple(float(value) for value in values)
        if not stamps or len(stamps) != len(levels):
            raise ValueError("times and values must be equally long and not empty")
        if not np.all(np.isfinite(stamps + levels)):
            raise ValueError("times and values must be finite")
        if np.any(np.diff(stamps) <= 0):
            raise ValueError("times must be strictly ascending")
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
    """

    electrode: Electrode
    waveform: float | Callable[[float], float]
    resistance: float = 0.0

    def __post_init__(self) -> None:
        require_nonnegative("resistance", self.resistance)

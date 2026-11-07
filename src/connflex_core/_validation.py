from collections.abc import Iterable

import numpy as np


def require_positive(name: str, value: float) -> None:
    if not value > 0:
        raise ValueError(f"{name} must be positive, got {value}")


def require_nonnegative(name: str, value: float) -> None:
    if not value >= 0:
        raise ValueError(f"{name} must be zero or positive, got {value}")


def require_at_least(name: str, value: float, bound: float) -> None:
    if not value >= bound:
        raise ValueError(f"{name} must be at least {bound}, got {value}")


def require_ascending(name: str, values: Iterable[float], least: int) -> tuple[float, ...]:
    """``values`` as floats: finite, strictly ascending, and at least ``least`` of them."""
    stamps = tuple(float(value) for value in values)
    if len(stamps) < least:
        raise ValueError(f"{name} must hold at least {least} value{'s' if least != 1 else ''}")
    if not np.all(np.isfinite(stamps)):
        raise ValueError(f"{name} must be finite")
    if np.any(np.diff(stamps) <= 0):
        raise ValueError(f"{name} must be strictly ascending")
    return stamps

def require_positive(name: str, value: float) -> None:
    if not value > 0:
        raise ValueError(f"{name} must be positive, got {value}")


def require_nonnegative(name: str, value: float) -> None:
    if not value >= 0:
        raise ValueError(f"{name} must be zero or positive, got {value}")

"""Values carried with their gradient and Hessian through arithmetic."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

Array = NDArray[np.float64]


@dataclass(frozen=True, eq=False)
class Jet:
    """Arrays of values with their first and second derivatives in d parameters, (103).

    Attributes:
        value: m values.
        gradient: m by d first derivatives.
        hessian: m by d by d second derivatives.
    """

    value: Array
    gradient: Array
    hessian: Array

    __array_ufunc__ = None  # an array on the left of an operator defers to the methods here

    @classmethod
    def constant(cls, value: ArrayLike, parameters: int) -> Jet:
        """Values that do not depend on any parameter."""
        value = np.asarray(value, dtype=float)
        return cls(
            value,
            np.zeros((*value.shape, parameters)),
            np.zeros((*value.shape, parameters, parameters)),
        )

    def __add__(self, other: Jet | float | Array) -> Jet:
        if isinstance(other, Jet):
            return Jet(
                self.value + other.value,
                self.gradient + other.gradient,
                self.hessian + other.hessian,
            )
        return Jet(self.value + other, self.gradient, self.hessian)

    __radd__ = __add__

    def __neg__(self) -> Jet:
        return Jet(-self.value, -self.gradient, -self.hessian)

    def __sub__(self, other: Jet | float | Array) -> Jet:
        return self + (-other)

    def __rsub__(self, other: float | Array) -> Jet:
        return -self + other

    def __mul__(self, other: Jet | float | Array) -> Jet:
        if isinstance(other, Jet):
            cross = self.gradient[:, :, None] * other.gradient[:, None, :]
            return Jet(
                self.value * other.value,
                self.value[:, None] * other.gradient + other.value[:, None] * self.gradient,
                self.value[:, None, None] * other.hessian
                + other.value[:, None, None] * self.hessian
                + cross
                + cross.transpose(0, 2, 1),
            )
        factor = np.asarray(other, dtype=float)
        return Jet(
            self.value * factor,
            self.gradient * factor[..., None],
            self.hessian * factor[..., None, None],
        )

    __rmul__ = __mul__

    def __truediv__(self, other: Jet | float | Array) -> Jet:
        if isinstance(other, Jet):
            return self * other.through(lambda x: 1 / x, lambda x: -1 / x**2, lambda x: 2 / x**3)
        return self * (1 / np.asarray(other, dtype=float))

    def through(
        self,
        function: Callable[[Array], Array],
        first: Callable[[Array], Array],
        second: Callable[[Array], Array],
    ) -> Jet:
        """A function of these values, given with its first and second derivative."""
        slope, bend = first(self.value), second(self.value)
        outer = self.gradient[:, :, None] * self.gradient[:, None, :]
        return Jet(
            function(self.value),
            slope[:, None] * self.gradient,
            slope[:, None, None] * self.hessian + bend[:, None, None] * outer,
        )

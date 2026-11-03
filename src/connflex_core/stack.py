"""Wafer layer stack."""

from __future__ import annotations

from dataclasses import dataclass

from connflex_core._validation import require_positive
from connflex_core.materials import Material


@dataclass(frozen=True)
class Layer:
    """One layer of the stack.

    Attributes:
        thickness: Layer thickness in meters.
        material: Material of the layer.
    """

    thickness: float
    material: Material

    def __post_init__(self) -> None:
        require_positive("thickness", self.thickness)


@dataclass(frozen=True)
class Stack:
    """SOI layer stack, listed from the top down.

    Attributes:
        device: Device layer. Its material also sets the beam stiffness.
        box: Buried oxide layer.
        handle: Handle layer.
    """

    device: Layer
    box: Layer
    handle: Layer

    @property
    def layers(self) -> tuple[Layer, Layer, Layer]:
        """Layers from the top down."""
        return (self.device, self.box, self.handle)

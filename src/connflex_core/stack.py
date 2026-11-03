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
        material: Layer material.
    """

    thickness: float
    material: Material

    def __post_init__(self) -> None:
        require_positive("thickness", self.thickness)


@dataclass(frozen=True)
class Stack:
    """SOI layer stack, listed from top down.

    Attributes:
        device: Device layer. Its material also sets beam stiffness.
        box: Buried oxide layer.
        handle: Handle layer.
    """

    device: Layer
    box: Layer
    handle: Layer

    @property
    def layers(self) -> tuple[Layer, Layer, Layer]:
        """Layers from top down."""
        return (self.device, self.box, self.handle)

    @property
    def z_extents(self) -> tuple[tuple[float, float], ...]:
        """Top and bottom height of each layer, with the device layer's top at z = 0."""
        extents = []
        top = 0.0
        for layer in self.layers:
            extents.append((top, top - layer.thickness))
            top -= layer.thickness
        return tuple(extents)

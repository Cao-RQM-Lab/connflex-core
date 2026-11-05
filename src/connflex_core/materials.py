"""Material properties."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from connflex_core._validation import require_nonnegative, require_positive

TORSION_MODELS = ("rectangle",)


@dataclass(frozen=True)
class Material:
    """Isotropic material of one stack layer.

    Attributes:
        density: Mass density in kilograms per cubic meter.
        youngs_modulus: Young's modulus in pascals.
        poisson_ratio: Poisson's ratio, used to derive shear modulus when it is not given.
        shear_modulus: Shear modulus in pascals. Derived as E / (2 (1 + nu)) when omitted.
        torsion_model: Torsion constant model for beams. "rectangle" is the Saint-Venant
            constant of a solid rectangle.
        torsion_beta: Short-beam torsional stiffening coefficient. Zero disables the correction.
        resistivity: Electrical resistivity in ohm meters. None leaves resistance of beam
            links unknown. Only the device layer's value is used.
    """

    density: float
    youngs_modulus: float
    poisson_ratio: float
    shear_modulus: float | None = None
    torsion_model: str = "rectangle"
    torsion_beta: float = 0.63
    resistivity: float | None = None

    def __post_init__(self) -> None:
        require_positive("density", self.density)
        require_positive("youngs_modulus", self.youngs_modulus)
        if not 0.0 <= self.poisson_ratio < 0.5:
            raise ValueError(f"poisson_ratio must be in [0, 0.5), got {self.poisson_ratio}")
        if self.shear_modulus is None:
            derived = self.youngs_modulus / (2.0 * (1.0 + self.poisson_ratio))
            object.__setattr__(self, "shear_modulus", derived)
        else:
            require_positive("shear_modulus", self.shear_modulus)
        if self.torsion_model not in TORSION_MODELS:
            raise ValueError(
                f"torsion_model must be one of {TORSION_MODELS}, got {self.torsion_model!r}"
            )
        require_nonnegative("torsion_beta", self.torsion_beta)
        if self.resistivity is not None:
            require_positive("resistivity", self.resistivity)

    @classmethod
    def silicon(cls, **overrides: Any) -> Material:
        """Single-crystal silicon: 2330 kg/m^3, E = 169 GPa along <110>, nu = 0.28."""
        return cls(
            **{"density": 2330.0, "youngs_modulus": 169e9, "poisson_ratio": 0.28, **overrides}
        )

    @classmethod
    def oxide(cls, **overrides: Any) -> Material:
        """Thermal silicon dioxide: 2200 kg/m^3, E = 70 GPa, nu = 0.17."""
        return cls(
            **{"density": 2200.0, "youngs_modulus": 70e9, "poisson_ratio": 0.17, **overrides}
        )

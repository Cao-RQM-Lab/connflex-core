"""Pair integrals and derivatives of the field kernel: compiled when numba is installed,
numpy otherwise."""

from __future__ import annotations

try:
    from connflex_core.electrostatics.compiled import pair_derivatives, pair_integrals

    COMPILED = True
except ImportError:
    from connflex_core.electrostatics.kernel import pair_derivatives, pair_integrals

    COMPILED = False

__all__ = ["COMPILED", "pair_derivatives", "pair_integrals"]

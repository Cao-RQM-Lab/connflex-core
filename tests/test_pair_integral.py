"""The pair integral against quadrature, with or without numba."""

import numpy as np
import pytest
import scipy.integrate

from connflex_core.electrostatics import kernel, pairs
from connflex_core.electrostatics.kernel import SERIES_REACH


def exact(a0: complex, a1: complex, b0: complex, b1: complex) -> float:
    """Double integral of ln|r - r'| over two elements by adaptive quadrature."""
    la, lb = abs(a1 - a0), abs(b1 - b0)

    def integrand(t: float, s: float) -> float:
        return float(np.log(abs((a0 + s * (a1 - a0)) - (b0 + t * (b1 - b0)))))

    value, _ = scipy.integrate.dblquad(integrand, 0.0, 1.0, 0.0, 1.0, epsabs=1e-14, epsrel=1e-13)
    return value * la * lb


@pytest.mark.parametrize("integrals", [kernel.pair_integrals, pairs.pair_integrals])
def test_pair_integral_matches_quadrature_on_both_sides_of_the_switch(integrals) -> None:  # type: ignore[no-untyped-def]
    a0, a1 = 0.0 + 0.0j, 1.0 + 0.0j
    threshold = (1.0 + 0.6) / (2 * SERIES_REACH)
    cases = [
        (0.3 + 0.4j, 0.9 + 0.7j),  # near, closed form
        (3.0 + 2.0j, 3.0 + 2.6j),  # within the switch, closed form
        (threshold + 0.5 - 0.01 + 1e-9j, threshold + 0.5 + 0.59 + 1e-9j),  # just inside
        (threshold + 0.5 + 0.01 + 1e-9j, threshold + 0.5 + 0.61 + 1e-9j),  # just outside: series
        (40.0 + 25.0j, 40.0 + 25.6j),  # far, series
    ]
    for b0, b1 in cases:
        found = float(integrals(*(np.array([z]) for z in (a0, a1, b0, b1)))[0])
        assert found == pytest.approx(
            exact(a0, a1, b0, b1), abs=1e-12 * abs(a1 - a0) * abs(b1 - b0)
        )


@pytest.mark.parametrize("integrals", [kernel.pair_integrals, pairs.pair_integrals])
def test_pair_integral_of_an_element_with_itself_and_its_neighbor(integrals) -> None:  # type: ignore[no-untyped-def]
    a0, a1 = 0.2 + 0.1j, 1.4 + 0.1j
    length = abs(a1 - a0)
    own = float(integrals(*(np.array([z]) for z in (a0, a1, a0, a1)))[0])
    assert own == pytest.approx(length**2 * (np.log(length) - 1.5), rel=1e-12)  # equation (68)
    b1 = 1.4 + 1.0j
    shared = float(integrals(*(np.array([z]) for z in (a0, a1, a1, b1)))[0])
    assert shared == pytest.approx(exact(a0, a1, a1, b1), abs=1e-9 * length * abs(b1 - a1))

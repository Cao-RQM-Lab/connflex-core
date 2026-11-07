import numpy as np
import pytest
from shapely import MultiPolygon, box

from connflex_core.electrostatics import kernel, pairs
from connflex_core.electrostatics.outlines import mesh_outlines

pytest.importorskip("numba")
from connflex_core.electrostatics import compiled


@pytest.fixture(scope="module")
def elements() -> tuple[np.ndarray, ...]:
    """Elements of two squares 0.2 apart; their pairs are identical, touching, near, and far."""
    squares = [MultiPolygon([box(0, 0, 1, 1)]), MultiPolygon([box(1.2, 0.3, 2.2, 1.3)])]
    mesh = mesh_outlines(squares, 0.1, 4.0)
    start, end = (points @ (1, 1j) for points in (mesh.start, mesh.end))
    pivots = np.array([0.4 + 0.5j, 1.7 + 0.8j])[mesh.conductor]
    return start, end, pivots


def test_compiled_kernels_are_chosen_when_numba_is_installed() -> None:
    assert pairs.COMPILED
    assert pairs.pair_integrals is compiled.pair_integrals
    assert pairs.pair_derivatives is compiled.pair_derivatives


def test_compiled_integrals_match_numpy_on_every_kind_of_pair(
    elements: tuple[np.ndarray, ...],
) -> None:
    start, end, _ = elements
    expected = kernel.pair_integrals(start[:, None], end[:, None], start, end)
    found = compiled.pair_integrals(start[:, None], end[:, None], start, end)
    assert found.shape == (len(start), len(start))
    assert found == pytest.approx(expected, rel=1e-12, abs=1e-14 * np.max(np.abs(expected)))


def test_compiled_derivatives_match_numpy_on_every_kind_of_pair(
    elements: tuple[np.ndarray, ...],
) -> None:
    start, end, pivots = elements
    apart = start[:, None] != start
    ends = (start[:, None], end[:, None], start, end)
    expected = kernel.pair_derivatives(*ends, pivots[:, None], pivots)
    found = compiled.pair_derivatives(*ends, pivots[:, None], pivots)
    assert found[0].shape == (6, len(start), len(start))
    assert found[1].shape == (6, 6, len(start), len(start))
    for mine, theirs in zip(found, expected, strict=True):
        mine, theirs = mine[..., apart], theirs[..., apart]
        assert mine == pytest.approx(theirs, rel=1e-12, abs=1e-14 * np.max(np.abs(theirs)))


def test_scalars_broadcast_against_arrays(elements: tuple[np.ndarray, ...]) -> None:
    start, end, pivots = elements
    one = compiled.pair_integrals(start[3], end[3], start, end)
    assert one.shape == (len(start),)
    assert one == pytest.approx(kernel.pair_integrals(start[3], end[3], start, end), rel=1e-12)
    slope, bend = compiled.pair_derivatives(
        start[0], end[0], start[5], end[5], pivots[0], pivots[5]
    )
    assert slope.shape == (6,)
    assert bend.shape == (6, 6)

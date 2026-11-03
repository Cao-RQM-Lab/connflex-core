from collections.abc import Callable
from pathlib import Path

import numpy as np
import pytest
from shapely import box

from connflex_core import ExtractionParams, Layout, SolveOptions, Stack, Structure, extract
from connflex_core.mechanics.stiffness import effective_torsion_constant

UM = 1e-6
PAD_A = box(0, 0, 100, 100)
PAD_B = box(300, 0, 400, 100)
BEAM = box(100, 45, 300, 55)
HANDLE_A = box(-50, -50, 150, 150)
Params = Callable[..., ExtractionParams]
LayoutUm = Callable[..., Layout]


def single_beam(stack: Stack, params: Params, layout_um: LayoutUm) -> Structure:
    return extract(layout_um([PAD_A, BEAM, PAD_B], [HANDLE_A]), stack, params())


def closed_form_axial_and_torsion(structure: Structure) -> tuple[float, float]:
    (beam,) = structure.beams
    body = structure.bodies[1]
    material = structure.material
    k_axial = material.youngs_modulus * beam.width * beam.thickness / beam.length
    j = effective_torsion_constant(beam.width, beam.thickness, beam.length, material.torsion_beta)
    assert material.shear_modulus is not None
    k_torsion = material.shear_modulus * j / beam.length
    f_axial = np.sqrt(k_axial / body.inertial.mass) / (2 * np.pi)
    f_torsion = np.sqrt(k_torsion / body.inertial.inertia[0, 0]) / (2 * np.pi)
    return f_axial, f_torsion


def test_single_pad_on_a_beam_has_the_closed_form_axial_and_torsion_modes(
    stack: Stack, params: Params, layout_um: LayoutUm
) -> None:
    structure = single_beam(stack, params, layout_um)
    result = structure.modal()
    assert result.frequencies_hz.shape == (6,)
    assert result.modes.shape == (6, 6)
    assert result.dofs == tuple((1, name) for name in ("ux", "uy", "uz", "tx", "ty", "tz"))
    assert np.all(np.diff(result.frequencies_hz) >= 0)
    for expected in closed_form_axial_and_torsion(structure):
        assert np.min(np.abs(result.frequencies_hz / expected - 1)) < 1e-9
    translations = np.abs(result.modes[:3, 0])
    assert result.dofs[int(np.argmax(translations))] == (1, "uy")


def test_parallel_beams_add_their_stiffness(
    stack: Stack, params: Params, layout_um: LayoutUm
) -> None:
    one = single_beam(stack, params, layout_um)
    two = extract(
        layout_um([PAD_A, box(100, 30, 300, 40), box(100, 60, 300, 70), PAD_B], [HANDLE_A]),
        stack,
        params(),
    )
    assert len(two.beams) == 2
    f_axial_one, _ = closed_form_axial_and_torsion(one)
    assert np.min(np.abs(two.modal().frequencies_hz / (np.sqrt(2) * f_axial_one) - 1)) < 1e-9


def test_static_axial_load_deflects_by_the_axial_compliance(
    stack: Stack, params: Params, layout_um: LayoutUm
) -> None:
    structure = single_beam(stack, params, layout_um)
    (beam,) = structure.beams
    force = 1e-6
    displacement = structure.static({1: [force, 0, 0, 0, 0, 0]})
    k_axial = structure.material.youngs_modulus * beam.width * beam.thickness / beam.length
    stretch = force / k_axial
    assert displacement[0] == pytest.approx(stretch, rel=1e-9)
    assert np.all(np.abs(displacement[1:3]) < 1e-9 * stretch)
    assert np.all(np.abs(displacement[3:]) * beam.length < 1e-9 * stretch)


def test_sparse_and_dense_solvers_agree_on_the_lowest_modes(
    stack: Stack, params: Params, layout_um: LayoutUm
) -> None:
    chain = [PAD_A, BEAM, PAD_B, box(400, 45, 600, 55), box(600, 0, 700, 100)]
    structure = extract(layout_um(chain, [HANDLE_A]), stack, params())
    dense = structure.modal(SolveOptions(n_modes=3))
    sparse = structure.modal(SolveOptions(n_modes=3, method="sparse"))
    assert dense.frequencies_hz.shape == sparse.frequencies_hz.shape == (3,)
    assert sparse.frequencies_hz == pytest.approx(dense.frequencies_hz, rel=1e-6)
    assert dense.mass.shape == dense.stiffness.shape == (12, 12)


def test_result_round_trips_through_npz(
    stack: Stack, params: Params, layout_um: LayoutUm, tmp_path: Path
) -> None:
    result = single_beam(stack, params, layout_um).modal()
    path = tmp_path / "modes.npz"
    result.to_npz(path)
    with np.load(path) as saved:
        assert saved["frequencies_hz"] == pytest.approx(result.frequencies_hz)
        assert saved["modes"].shape == (6, 6)
        assert list(saved["dof_name"]) == ["ux", "uy", "uz", "tx", "ty", "tz"]

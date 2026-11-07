import pytest

from connflex_core import ElectrostaticParams, ExtractionParams, SolveOptions


def params(**overrides: object) -> ExtractionParams:
    values: dict[str, object] = {
        "release_hole_max_area": 100e-12,
        "release_distance": 5e-6,
        "body_min_width": 10e-6,
    }
    values.update(overrides)
    return ExtractionParams(**values)  # type: ignore[arg-type]


def test_dependent_defaults_follow_their_source_fields() -> None:
    p = params(heal_distance=20e-9)
    assert p.simplify_tolerance == 20e-9
    assert p.beam_overlap_margin == 20e-9
    assert p.beam_width_max == 10e-6


def test_explicit_values_override_dependent_defaults() -> None:
    p = params(simplify_tolerance=5e-9, beam_overlap_margin=1e-6, beam_width_max=8e-6)
    assert (p.simplify_tolerance, p.beam_overlap_margin, p.beam_width_max) == (5e-9, 1e-6, 8e-6)


@pytest.mark.parametrize(
    ("overrides", "field"),
    [
        ({"release_distance": 0.0}, "release_distance"),
        ({"release_hole_max_area": -1.0}, "release_hole_max_area"),
        ({"beam_overlap_margin": 5e-6}, "beam_overlap_margin"),
        ({"beam_width_min": 10e-6}, "beam_width_min"),
        ({"buffer_resolution": 0}, "buffer_resolution"),
        ({"precision_grid": 0.0}, "precision_grid"),
        ({"anchor": (1.0,)}, "anchor"),
        ({"via_resistance": -1.0}, "via_resistance"),
        ({"heal_distance": 0.0}, "heal_distance"),
        ({"simplify_tolerance": 0.0}, "simplify_tolerance"),
        ({"body_min_width": 0.0}, "body_min_width"),
        ({"beam_width_max": 0.0}, "beam_width_max"),
        ({"anchor": -1}, "anchor"),
    ],
)
def test_invalid_parameters_name_the_field(overrides: dict[str, object], field: str) -> None:
    with pytest.raises(ValueError, match=field):
        params(**overrides)


def test_anchor_accepts_index_or_point() -> None:
    assert params(anchor=2).anchor == 2
    assert params(anchor=(1e-6, 2e-6)).anchor == (1e-6, 2e-6)


def test_solve_options_reject_unknown_method_and_no_modes() -> None:
    with pytest.raises(ValueError, match="method"):
        SolveOptions(method="magic")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="n_modes"):
        SolveOptions(n_modes=0)


@pytest.mark.parametrize(
    ("overrides", "field"),
    [
        ({"ambient_permittivity": 0.5}, "ambient_permittivity"),
        ({"element_size": 0.0}, "element_size"),
        ({"corner_grading": 0.5}, "corner_grading"),
        ({"max_elements": 0}, "max_elements"),
        ({"interaction_distance": 0.0}, "interaction_distance"),
        ({"window_elements": 0}, "window_elements"),
        ({"window_tolerance": 0.0}, "window_tolerance"),
    ],
)
def test_invalid_electrostatic_parameters_name_the_field(
    overrides: dict[str, float], field: str
) -> None:
    with pytest.raises(ValueError, match=field):
        ElectrostaticParams(**overrides)  # type: ignore[arg-type]

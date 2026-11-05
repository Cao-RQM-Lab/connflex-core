import numpy as np
import pytest

from connflex_core import Layer, Material, Stack


def test_silicon_preset_derives_shear_modulus() -> None:
    silicon = Material.silicon()
    assert silicon.shear_modulus == pytest.approx(169e9 / (2 * (1 + 0.28)))


def test_preset_accepts_overrides() -> None:
    assert Material.silicon(youngs_modulus=130e9).youngs_modulus == 130e9
    assert Material.oxide(torsion_beta=0.0).torsion_beta == 0.0


def test_explicit_shear_modulus_is_kept() -> None:
    assert Material.silicon(shear_modulus=50e9).shear_modulus == 50e9


@pytest.mark.parametrize(
    ("overrides", "field"),
    [
        ({"density": -1.0}, "density"),
        ({"youngs_modulus": 0.0}, "youngs_modulus"),
        ({"poisson_ratio": 0.5}, "poisson_ratio"),
        ({"torsion_model": "circle"}, "torsion_model"),
        ({"torsion_beta": -0.1}, "torsion_beta"),
        ({"resistivity": 0.0}, "resistivity"),
    ],
)
def test_invalid_material_values_name_the_field(overrides: dict[str, object], field: str) -> None:
    with pytest.raises(ValueError, match=field):
        Material.silicon(**overrides)


def test_layer_requires_positive_thickness() -> None:
    with pytest.raises(ValueError, match="thickness"):
        Layer(thickness=0.0, material=Material.silicon())


def test_stack_orders_layers_from_the_top_down() -> None:
    stack = Stack(
        device=Layer(thickness=25e-6, material=Material.silicon()),
        box=Layer(thickness=2e-6, material=Material.oxide()),
        handle=Layer(thickness=400e-6, material=Material.silicon()),
    )
    assert [layer.thickness for layer in stack.layers] == [25e-6, 2e-6, 400e-6]
    expected = np.array([(0.0, -25e-6), (-25e-6, -27e-6), (-27e-6, -427e-6)])
    assert np.asarray(stack.z_extents) == pytest.approx(expected)

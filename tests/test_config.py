import re
from pathlib import Path

import pytest

from connflex_core import load_config

PAGE = Path(__file__).parents[1] / "docs" / "parameters.md"


def documented_example() -> str:
    return re.search(r"```toml\n(.*?)```", PAGE.read_text(), flags=re.S).group(1)  # type: ignore[union-attr]


def test_documented_example_loads(tmp_path: Path) -> None:
    path = tmp_path / "process.toml"
    path.write_text(documented_example())
    stack, params, options = load_config(path)
    assert stack.device.thickness == 25e-6
    assert stack.box.material.density == 2200.0
    assert stack.handle.material.youngs_modulus == 130e9
    assert params.body_min_width == 10e-6
    assert options.n_modes == 20


def test_unknown_key_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "process.toml"
    path.write_text(documented_example().replace("release_distance", "release_distnace"))
    with pytest.raises(ValueError, match="release_distnace"):
        load_config(path)


def test_unknown_material_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "process.toml"
    path.write_text(documented_example().replace('material = "oxide"', 'material = "glass"'))
    with pytest.raises(ValueError, match="glass"):
        load_config(path)


def test_new_material_needs_its_required_properties(tmp_path: Path) -> None:
    path = tmp_path / "process.toml"
    text = documented_example().replace('material = "oxide"', 'material = "nitride"')
    path.write_text(text + "\n[materials.nitride]\ndensity = 3100\n")
    with pytest.raises(ValueError, match="youngs_modulus"):
        load_config(path)

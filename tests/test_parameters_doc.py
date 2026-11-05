"""Parameters page and configuration classes must list the same fields."""

import dataclasses
import re
from pathlib import Path

import pytest

from connflex_core import (
    ElectrostaticParams,
    ExtractionParams,
    Layer,
    Material,
    SolveOptions,
    StaticOptions,
    Wire,
)

PAGE = Path(__file__).parents[1] / "docs" / "parameters.md"


def documented_rows(heading_contains: str) -> dict[str, set[str]]:
    """Return {class prefix or '': field names} for the table under the matching heading."""
    sections = re.split(r"^## ", PAGE.read_text(), flags=re.M)
    section = next(s for s in sections if heading_contains in s.splitlines()[0])
    rows: dict[str, set[str]] = {}
    for prefix, name in re.findall(r"^\| `(?:(\w+)\.)?(\w+)` \|", section, flags=re.M):
        rows.setdefault(prefix, set()).add(name)
    return rows


def field_names(cls: type) -> set[str]:
    return {f.name for f in dataclasses.fields(cls)}


def test_extraction_params_match_the_page() -> None:
    assert documented_rows("`ExtractionParams`")[""] == field_names(ExtractionParams)


def test_solve_options_match_the_page() -> None:
    assert documented_rows("`SolveOptions`")[""] == field_names(SolveOptions)


def test_electrostatic_params_match_the_page() -> None:
    assert documented_rows("`ElectrostaticParams`")[""] == field_names(ElectrostaticParams)


def test_static_options_match_the_page() -> None:
    assert documented_rows("`StaticOptions`")[""] == field_names(StaticOptions)


def test_wire_matches_the_page() -> None:
    assert documented_rows("`Wire`")[""] == field_names(Wire)


def test_stack_and_material_match_the_page() -> None:
    rows = documented_rows("`Material`")
    assert rows["Material"] == field_names(Material)
    assert rows["Layer"] == field_names(Layer) - {"material"}


@pytest.mark.parametrize(
    "cls",
    [ExtractionParams, SolveOptions, ElectrostaticParams, StaticOptions, Wire, Layer, Material],
)
def test_every_field_is_described_in_the_class_docstring(cls: type) -> None:
    doc = cls.__doc__ or ""
    attributes = doc.split("Attributes:")[1]
    described = set(re.findall(r"^\s+(\w+):", attributes, flags=re.M))
    assert described == field_names(cls)

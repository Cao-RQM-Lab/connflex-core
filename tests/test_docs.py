"""The documentation runs and stays consistent: the quick start and the usage example
execute, theory equations form one numbered sequence, and every citation of one exists."""

import re
import textwrap
from pathlib import Path

import numpy as np
import pytest

DOCS = Path(__file__).parents[1] / "docs"


def python_blocks(text: str) -> list[str]:
    """Fenced Python blocks of a page, dedented, in order."""
    found = re.findall(r"^( *)```python\n(.*?)^\1```", text, flags=re.M | re.S)
    return [textwrap.dedent(body) for _, body in found]


def block_after(text: str, marker: str) -> str:
    return python_blocks(text[text.index(marker) :])[0]


def test_quick_start_runs_through_all_three_guides() -> None:
    namespace: dict[str, object] = {}
    for block in python_blocks((DOCS / "quickstart.md").read_text()):
        exec(block, namespace)
    structure = namespace["structure"]
    assert len(structure.bodies) == 2  # type: ignore[attr-defined]
    assert len(structure.beams) == 4  # type: ignore[attr-defined]
    assert structure.bodies[1].box.is_empty  # type: ignore[attr-defined]
    assert namespace["modes"].frequencies_hz.shape == (6,)  # type: ignore[attr-defined]
    assert float(namespace["gap_capacitance"]) > 0  # type: ignore[arg-type]
    assert float(namespace["force"][namespace["uz"]]) < 0  # type: ignore[index]
    assert float(namespace["softening"]) < 0  # type: ignore[arg-type]
    assert namespace["walls"].method[0, 1] in ("dense", "windowed")  # type: ignore[attr-defined]
    bias, run, uz = namespace["bias"], namespace["run"], namespace["uz"]
    assert bias.stable  # type: ignore[attr-defined]
    assert 0 < bias.margin < 1  # type: ignore[attr-defined]
    assert bias.displacement[uz] < 0  # type: ignore[attr-defined]
    assert namespace["f_biased"] < namespace["f_vertical"]  # type: ignore[operator]
    assert np.isfinite(namespace["admittance"]).all()
    assert "plate" in namespace["netlist"]  # type: ignore[operator]
    assert run.displacement[-1, uz] == pytest.approx(bias.displacement[uz], rel=1e-2)  # type: ignore[attr-defined]
    assert abs(run.energy.imbalance[-1]) < 1e-3 * abs(run.energy.supplied[-1])  # type: ignore[attr-defined]


def test_usage_example_extracts_a_movable_body_and_its_modes() -> None:
    text = (DOCS / "usage.md").read_text()
    namespace: dict[str, object] = {}
    for marker in ('=== "shapely"', "## Describe the process", "## Choose parameters"):
        exec(block_after(text, marker), namespace)
    cf = namespace["cf"]
    structure = cf.extract(namespace["layout"], namespace["stack"], namespace["params"])  # type: ignore[attr-defined]
    assert len(structure.beams) == 1
    assert [body.anchor for body in structure.bodies] == [True, False]
    assert not structure.rejected
    result = structure.modal(namespace["options"])
    assert result.frequencies_hz.shape == (6,)
    assert result.frequencies_hz[0] > 0


def theory_pages_in_order() -> list[Path]:
    nav = (DOCS.parent / "mkdocs.yml").read_text()
    return [DOCS / name for name in re.findall(r"(theory/[\w-]+\.md)", nav)]


def test_theory_equations_form_one_sequence_across_the_pages() -> None:
    pages = theory_pages_in_order()
    assert {page.name for page in pages} == {page.name for page in (DOCS / "theory").glob("*.md")}
    tags: list[str] = []
    for page in pages:
        tags += re.findall(r"\\tag\{(\d+[a-c]?)\}", page.read_text())
    assert len(tags) == len(set(tags))
    numbers = [int(re.match(r"\d+", tag).group()) for tag in tags]  # type: ignore[union-attr]
    assert numbers == sorted(numbers)
    assert set(numbers) == set(range(1, max(numbers) + 1))


def test_every_cited_equation_exists() -> None:
    tags: set[str] = set()
    for page in theory_pages_in_order():
        tags |= set(re.findall(r"\\tag\{(\d+[a-c]?)\}", page.read_text()))
    for page in [*DOCS.glob("*.md"), *(DOCS / "theory").glob("*.md")]:
        cited = set(re.findall(r"\((\d{1,3}[a-c]?)\)", page.read_text()))
        numbered = {tag for tag in cited if 1 <= int(re.match(r"\d+", tag).group()) <= 200}  # type: ignore[union-attr]
        missing = numbered - tags
        assert not missing, (page.name, sorted(missing))

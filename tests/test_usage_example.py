"""The introductory example on the usage page runs through to modes."""

import re
import textwrap
from pathlib import Path

PAGE = Path(__file__).parents[1] / "docs" / "usage.md"


def block_after(text: str, marker: str) -> str:
    """First fenced Python block after ``marker``, dedented."""
    tail = text[text.index(marker) :]
    match = re.search(r"^( *)```python\n(.*?)^\1```", tail, flags=re.M | re.S)
    assert match is not None, marker
    return textwrap.dedent(match.group(2))


def test_shapely_example_extracts_a_movable_body_and_its_modes() -> None:
    text = PAGE.read_text()
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

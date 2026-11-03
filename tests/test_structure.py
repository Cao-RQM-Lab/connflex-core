import json
from collections.abc import Callable

from shapely import box

from connflex_core import ExtractionParams, Layout, Stack, extract


def test_graph_and_tree_describe_the_structure(
    stack: Stack,
    params: Callable[..., ExtractionParams],
    layout_um: Callable[..., Layout],
) -> None:
    device = [box(0, 0, 100, 100), box(100, 45, 300, 55), box(300, 0, 400, 100)]
    structure = extract(layout_um(device, [box(-50, -50, 150, 150)]), stack, params())
    graph = structure.graph
    assert sorted(graph.nodes) == [0, 1]
    assert graph.number_of_edges() == 1
    tree = json.loads(json.dumps(structure.to_dict()))
    assert [body["anchor"] for body in tree["bodies"]] == [True, False]
    assert tree["beams"][0]["bodies"] == [0, 1]
    assert tree["bodies"][0]["device"]["type"] == "MultiPolygon"

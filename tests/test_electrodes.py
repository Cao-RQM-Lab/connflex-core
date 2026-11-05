from collections.abc import Callable

import networkx as nx
import pytest
from shapely import box

from connflex_core import ExtractionParams, Layout, Stack, Structure, extract

UM = 1e-6
PAD_A = box(0, 0, 100, 100)
PAD_B = box(300, 0, 400, 100)
BEAM = box(100, 45, 300, 55)
HANDLE_A = box(-150, -50, 150, 150)
HANDLE_AB = box(-50, -50, 450, 150)
Params = Callable[..., ExtractionParams]
LayoutUm = Callable[..., Layout]


def nets_joined_by_links(structure: Structure) -> list[set[int]]:
    graph: nx.Graph[int] = nx.Graph()
    graph.add_nodes_from(electrode.index for electrode in structure.electrodes)
    graph.add_edges_from(link.electrodes for link in structure.links)
    return sorted(
        ({structure.electrodes[i].net for i in part} for part in nx.connected_components(graph)),
        key=min,
    )


def test_bodies_joined_by_a_beam_are_electrodes_joined_by_a_link(
    stack: Stack, params: Params, layout_um: LayoutUm
) -> None:
    structure = extract(layout_um([PAD_A, BEAM, PAD_B], [HANDLE_A]), stack, params())
    pad_a, pad_b, handle = structure.electrodes
    assert [(e.index, e.layer, e.body) for e in structure.electrodes] == [
        (0, "device", 0),
        (1, "device", 1),
        (2, "handle", 0),
    ]
    assert pad_a.conductor == pad_b.conductor != handle.conductor
    assert pad_a.net == pad_b.net == structure.beams[0].net != handle.net
    assert pad_b.region.area == pytest.approx(PAD_B.area * UM * UM)
    (link,) = structure.links
    assert (link.index, link.kind, link.electrodes, link.source) == (0, "beam", (0, 1), 0)
    assert nets_joined_by_links(structure) == sorted([{pad_a.net}, {handle.net}], key=min)


def test_stubs_belong_to_the_electrode_of_the_polygon_they_touch(
    stack: Stack, params: Params, layout_um: LayoutUm
) -> None:
    pad_c = box(0, 300, 100, 400)
    finger_a, finger_c = box(100, 45, 200, 55), box(100, 345, 200, 355)
    handle = box(-50, -50, 150, 450)
    structure = extract(layout_um([PAD_A, finger_a, pad_c, finger_c], [handle]), stack, params())
    assert len(structure.bodies) == 1
    on_a = structure.electrode_at((150 * UM, 50 * UM), "device")
    on_c = structure.electrode_at((150 * UM, 350 * UM), "device")
    assert on_a.index == structure.electrode_at((50 * UM, 50 * UM), "device").index
    assert on_a.index != on_c.index
    assert on_a.body == on_c.body == 0
    assert on_a.net != on_c.net
    assert [len(e.region.geoms) for e in (on_a, on_c)] == [2, 2]
    assert on_a.region.area == pytest.approx((PAD_A.area + finger_a.area) * UM * UM)
    assert structure.electrode_at((0.0, 0.0), "handle").layer == "handle"
    assert structure.links == ()
    with pytest.raises(ValueError, match="no device electrode"):
        structure.electrode_at((250 * UM, 200 * UM), "device")


def test_candidate_touching_several_polygons_is_a_bridge_between_their_electrodes(
    stack: Stack, params: Params, layout_um: LayoutUm
) -> None:
    structure = extract(layout_um([PAD_A, BEAM, PAD_B], [HANDLE_AB]), stack, params())
    assert len(structure.bodies) == 1
    assert structure.beams == ()
    pad_a, pad_b, handle = structure.electrodes
    assert [len(e.region.geoms) for e in (pad_a, pad_b)] == [1, 1]
    (link,) = structure.links
    assert (link.kind, link.electrodes, link.source) == ("bridge", (0, 1), 0)
    assert structure.rejected[link.source].polygon.area == pytest.approx(BEAM.area * UM * UM)
    assert nets_joined_by_links(structure) == sorted([{pad_a.net}, {handle.net}], key=min)


def test_via_links_the_electrodes_it_overlaps(
    stack: Stack, params: Params, layout_um: LayoutUm
) -> None:
    via = box(40, 40, 60, 60)
    structure = extract(layout_um([PAD_A, PAD_B], [HANDLE_AB], vias=[via]), stack, params())
    pad_a = structure.electrode_at((50 * UM, 50 * UM), "device")
    pad_b = structure.electrode_at((350 * UM, 50 * UM), "device")
    handle = structure.electrode_at((200 * UM, 120 * UM), "handle")
    (link,) = structure.links
    assert (link.kind, link.source) == ("via", 0)
    assert link.electrodes == tuple(sorted((pad_a.index, handle.index)))
    assert nets_joined_by_links(structure) == sorted([{pad_a.net}, {pad_b.net}], key=min)


def test_region_touching_no_body_is_an_electrode_of_the_body_whose_oxide_carries_it(
    stack: Stack, params: Params, layout_um: LayoutUm
) -> None:
    trace = box(300, 110, 400, 125)
    island = box(250, -50, 450, 150)
    structure = extract(layout_um([PAD_A, BEAM, PAD_B, trace], [HANDLE_A, island]), stack, params())
    assert len(structure.bodies) == 2
    on_trace = structure.electrode_at((350 * UM, 117 * UM), "device")
    on_pad_b = structure.electrode_at((350 * UM, 50 * UM), "device")
    assert on_trace.index != on_pad_b.index
    assert on_trace.body == on_pad_b.body == 1
    assert on_trace.net != on_pad_b.net
    assert [e.layer for e in structure.electrodes] == ["device"] * 3 + ["handle"] * 2
    assert len(nets_joined_by_links(structure)) == structure.nets.count

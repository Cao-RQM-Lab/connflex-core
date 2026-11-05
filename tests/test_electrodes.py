from collections.abc import Callable

import networkx as nx
import pytest
from shapely import box

from connflex_core import (
    ExtractionParams,
    Layer,
    Layout,
    Material,
    Stack,
    Structure,
    Wire,
    extract,
)

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
    assert link.resistance is None
    assert nets_joined_by_links(structure) == sorted([{pad_a.net}, {handle.net}], key=min)


def test_beam_link_resistance_is_resistivity_times_length_over_cross_section(
    stack: Stack, params: Params, layout_um: LayoutUm
) -> None:
    doped = Layer(thickness=25 * UM, material=Material.silicon(resistivity=1e-4))
    conducting = Stack(device=doped, box=stack.box, handle=stack.handle)
    structure = extract(layout_um([PAD_A, BEAM, PAD_B], [HANDLE_A]), conducting, params())
    (link,) = structure.links
    length = 200.01 * UM  # overlap centroids sit half the 10 nm margin inside each pad
    width = BEAM.area * UM * UM / length
    assert link.resistance == pytest.approx(1e-4 * length / (width * 25 * UM), rel=1e-9)


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
    assert link.resistance == 0.0
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
    assert (link.kind, link.source, link.resistance) == ("via", 0, 0.0)
    assert link.electrodes == tuple(sorted((pad_a.index, handle.index)))
    assert nets_joined_by_links(structure) == sorted([{pad_a.net}, {pad_b.net}], key=min)


def test_via_resistance_applies_between_each_device_and_handle_electrode_it_joins(
    stack: Stack, params: Params, layout_um: LayoutUm
) -> None:
    across = box(50, 40, 350, 60)
    layout = layout_um([PAD_A, PAD_B], [HANDLE_AB], vias=[across])
    structure = extract(layout, stack, params(via_resistance=5.0))
    pad_a = structure.electrode_at((50 * UM, 50 * UM), "device")
    pad_b = structure.electrode_at((350 * UM, 50 * UM), "device")
    handle = structure.electrode_at((200 * UM, 120 * UM), "handle")
    assert [(link.kind, link.source, link.resistance) for link in structure.links] == [
        ("via", 0, 5.0)
    ] * 2
    assert {link.electrodes for link in structure.links} == {
        (pad_a.index, handle.index),
        (pad_b.index, handle.index),
    }
    assert len(nets_joined_by_links(structure)) == structure.nets.count == 1


def test_wire_off_the_die_joins_pads_that_silicon_keeps_apart(
    stack: Stack, params: Params, layout_um: LayoutUm
) -> None:
    """A bond between two pads, and a bond from one of them to the handle beneath."""
    bond = Wire((50 * UM, 50 * UM), (350 * UM, 50 * UM), resistance=0.2)
    ground = Wire((60 * UM, 60 * UM), (200 * UM, 120 * UM, "handle"))
    layout = layout_um([PAD_A, PAD_B], [HANDLE_AB])
    alone = extract(layout, stack, params())
    assert alone.nets.count == 3
    structure = extract(layout, stack, params(), wires=[bond, ground])
    pad_a = structure.electrode_at((50 * UM, 50 * UM), "device")
    pad_b = structure.electrode_at((350 * UM, 50 * UM), "device")
    handle = structure.electrode_at((200 * UM, 120 * UM), "handle")
    assert structure.nets.count == 1
    assert pad_a.net == pad_b.net == handle.net
    assert [(link.kind, link.source, link.resistance) for link in structure.links] == [
        ("wire", 0, 0.2),
        ("wire", 1, 0.0),
    ]
    assert [link.electrodes for link in structure.links] == [
        tuple(sorted((pad_a.index, pad_b.index))),
        tuple(sorted((pad_a.index, handle.index))),
    ]
    first, second = structure.nets.wires
    assert (first.wire, first.conductors) == (bond, (pad_a.conductor, pad_b.conductor))
    assert second.conductors == (pad_a.conductor, handle.conductor)
    assert structure.nets.graph.has_edge(pad_a.conductor, pad_b.conductor)
    assert nets_joined_by_links(structure) == [{pad_a.net}]
    assert structure.to_dict()["links"][0]["kind"] == "wire"


def test_wire_ends_must_land_on_an_electrode(
    stack: Stack, params: Params, layout_um: LayoutUm
) -> None:
    layout = layout_um([PAD_A, BEAM, PAD_B, box(0, 300, 100, 400)], [box(-50, -50, 450, 450)])
    pad = (50 * UM, 50 * UM)
    with pytest.raises(ValueError, match=r"wire 1: .* no device conductor"):
        extract(
            layout, stack, params(), wires=[Wire(pad, (50 * UM, 350 * UM)), Wire(pad, (0.0, -1e-3))]
        )
    with pytest.raises(ValueError, match=r"wire 0: .* beam or bridge"):
        extract(layout, stack, params(), wires=[Wire((50 * UM, 350 * UM), (200 * UM, 50 * UM))])
    with pytest.raises(ValueError, match=r"wire 0: both ends lie on one electrode"):
        extract(layout, stack, params(), wires=[Wire(pad, (60 * UM, 60 * UM))])
    with pytest.raises(ValueError, match="resistance"):
        Wire(pad, (350 * UM, 50 * UM), resistance=-1.0)
    with pytest.raises(ValueError, match="layer"):
        Wire(pad, (350 * UM, 50 * UM, "box"))  # type: ignore[arg-type]


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

from collections.abc import Callable

import pytest
from shapely import box

from connflex_core import ExtractionParams, Layout, Stack, extract

UM = 1e-6
PAD_A = box(0, 0, 100, 100)
PAD_B = box(300, 0, 400, 100)
BEAM = box(100, 45, 300, 55)
HANDLE_AB = box(-50, -50, 450, 150)
HANDLE_A = box(-50, -50, 150, 150)
Params = Callable[..., ExtractionParams]
LayoutUm = Callable[..., Layout]


def test_isolated_pads_and_the_handle_are_separate_nets(
    stack: Stack, params: Params, layout_um: LayoutUm
) -> None:
    nets = extract(layout_um([PAD_A, PAD_B], [HANDLE_AB]), stack, params()).nets
    assert nets.count == 3
    assert [conductor.layer for conductor in nets.conductors] == ["device", "device", "handle"]
    assert len({conductor.net for conductor in nets.conductors}) == 3
    assert nets.graph.number_of_edges() == 0


def test_vias_merge_the_conductors_they_overlap(
    stack: Stack, params: Params, layout_um: LayoutUm
) -> None:
    on_a = box(40, 40, 60, 60)
    nets = extract(layout_um([PAD_A, PAD_B], [HANDLE_AB], vias=[on_a]), stack, params()).nets
    assert nets.count == 2
    pad_a = nets.conductor_at((50 * UM, 50 * UM), "device")
    pad_b = nets.conductor_at((350 * UM, 50 * UM), "device")
    handle = nets.conductor_at((200 * UM, 120 * UM), "handle")
    assert pad_a.net == handle.net != pad_b.net
    assert nets.graph.has_edge(pad_a.index, handle.index)

    across = box(50, 40, 350, 60)
    nets = extract(layout_um([PAD_A, PAD_B], [HANDLE_AB], vias=[across]), stack, params()).nets
    assert nets.count == 1


def test_beams_and_bodies_carry_the_net_of_their_conductor(
    stack: Stack, params: Params, layout_um: LayoutUm
) -> None:
    pad_c = box(0, 300, 100, 400)
    handle = box(-50, -50, 150, 450)
    structure = extract(layout_um([PAD_A, BEAM, PAD_B, pad_c], [handle]), stack, params())
    (beam,) = structure.beams
    anchor, free = structure.bodies
    assert structure.nets.count == 3
    assert free.nets == (beam.net,)
    assert beam.net in anchor.nets
    assert len(set(anchor.nets)) == 2


def test_via_without_a_handle_conductor_is_ignored_with_a_warning(
    stack: Stack, params: Params, layout_um: LayoutUm
) -> None:
    stray = box(340, 40, 360, 60)
    with pytest.warns(UserWarning, match="via 0"):
        structure = extract(
            layout_um([PAD_A, BEAM, PAD_B], [HANDLE_A], vias=[stray]), stack, params()
        )
    assert structure.nets.count == 2
    assert structure.nets.graph.number_of_edges() == 0


def test_find_nets_takes_bare_polygons() -> None:
    from connflex_core import Wire, find_nets

    device = [box(0, 0, 100e-6, 100e-6), box(300e-6, 0, 400e-6, 100e-6)]
    handle = [box(-50e-6, -50e-6, 450e-6, 150e-6)]
    nets = find_nets(device, handle)
    assert nets.count == 3
    joined = find_nets(device, handle, wires=[Wire((50e-6, 50e-6), (350e-6, 50e-6))])
    assert joined.count == 2
    assert joined.conductors[0].net == joined.conductors[1].net
    through = find_nets(device, handle, vias=[box(40e-6, 40e-6, 60e-6, 60e-6)])
    assert through.count == 2
    assert through.conductors[0].net == through.conductors[2].net

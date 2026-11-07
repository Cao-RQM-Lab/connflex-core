"""Layout to structure."""

from __future__ import annotations

from collections.abc import Sequence

import shapely
from shapely import MultiPolygon, Polygon, STRtree

from connflex_core.electrical.electrodes import find_electrodes
from connflex_core.electrical.nets import Wire, find_nets
from connflex_core.geometry.beams import beam_candidates, find_beams
from connflex_core.geometry.bodies import layered_components, select_anchor
from connflex_core.geometry.heal import as_multipolygon, heal, opening, snap
from connflex_core.geometry.release import fill_release_holes, openings_within, surviving_box
from connflex_core.layout import Layout
from connflex_core.mechanics.inertia import mass_properties
from connflex_core.model import Body, Structure
from connflex_core.params import ExtractionParams
from connflex_core.progress import Progress, stage
from connflex_core.stack import Stack


def extract(
    layout: Layout,
    stack: Stack,
    params: ExtractionParams,
    *,
    wires: Sequence[Wire] = (),
    progress: Progress | None = None,
) -> Structure:
    """Extract rigid bodies, beams, and electrical nets from a layout.

    Release holes inside the device layer are filled, the layers are healed and simplified,
    and holes that healing closes are filled as well. The holes and the release layer give
    the surviving buried oxide; holes are disjoint by construction and the release layer is
    unioned by its reader, so the two are unioned only when both are present. An anchor handle polygon is selected and every polygon
    with no connection to it through the stack is dropped. Bodies are assembled by overlap
    through the opened device layer, oxide, and handle; each gets mass properties from the
    stack, with the moments of the release openings inside it subtracted from its device
    layer. Device regions removed by
    the opening are qualified as beams when they touch exactly two bodies, keep no surviving
    oxide underneath, and lie within the width bounds; a candidate touching exactly one body
    is a stub that adds mass to that body. Remaining device and handle polygons are conductors of the
    nets, joined through vias and through ``wires``, connections made off the die between
    pads that the layout leaves apart. Each body polygon with the stubs that touch it
    alone is an electrode, and beams, candidates that touch several body polygons, vias,
    and wires are links between electrodes.

    Stages heal, release, prune, nets, opening, bodies, candidates, beams, mass, and
    electrodes each log through the ``connflex_core`` logger and report to ``progress``.
    """
    tolerance = params.simplify_tolerance
    margin = params.beam_overlap_margin
    width_max = params.beam_width_max
    if tolerance is None or margin is None or width_max is None:
        raise ValueError("extraction parameters are not resolved")

    with stage("heal", progress) as current:
        device, holes = fill_release_holes(layout.device, params.release_hole_max_area)
        device = heal(device, params.heal_distance, tolerance)
        device, closed = fill_release_holes(device, params.release_hole_max_area)
        holes += closed
        handle = heal(layout.handle, params.heal_distance, tolerance)
        release = None
        if layout.release is not None:
            release = as_multipolygon(layout.release.simplify(tolerance, preserve_topology=True))
        if params.precision_grid is not None:
            device = snap(device, params.precision_grid)
            handle = snap(handle, params.precision_grid)
            release = None if release is None else snap(release, params.precision_grid)
            if holes:
                holes = list(snap(MultiPolygon(holes), params.precision_grid).geoms)
        current.count(device=len(device.geoms), handle=len(handle.geoms), holes=len(holes))

    with stage("release", progress) as current:
        if release is None:
            openings = MultiPolygon(holes)
        elif not holes:
            openings = release
        else:
            openings = as_multipolygon(shapely.union_all([*holes, *release.geoms]))
        box = surviving_box(
            device, handle, openings, params.release_distance, params.buffer_resolution
        )
        opening_parts = list(openings.geoms)
        opening_tree = STRtree(opening_parts)
        current.count(openings=len(opening_parts), box=len(box.geoms))

    handle_parts = _by_area(handle)
    anchor = select_anchor(handle_parts, params.anchor)
    device_parts, box_parts = list(device.geoms), list(box.geoms)
    layers = (device_parts, box_parts, handle_parts)
    with stage("prune", progress, total=sum(map(len, layers))) as current:
        reachable = layered_components(layers, seed=(2, anchor), advance=current.advance)[0]
        kept = [sorted(i for layer, i in reachable if layer == which) for which in range(3)]
        device_parts = [device_parts[i] for i in kept[0]]
        box_parts = [box_parts[i] for i in kept[1]]
        handle_parts = [handle_parts[i] for i in kept[2]]
        anchor = kept[2].index(anchor)
        current.count(dropped=sum(map(len, layers)) - len(reachable))
    device = MultiPolygon(device_parts)

    with stage("nets", progress) as current:
        vias = [] if layout.vias is None else list(layout.vias.geoms)
        nets = find_nets(device_parts, handle_parts, vias, wires)
        current.count(nets=nets.count, vias=len(vias), wires=len(wires))
    conductors = STRtree(device_parts)

    def net_of(polygon: Polygon) -> int:
        return nets.conductors[int(conductors.nearest(polygon.representative_point()))].net

    with stage("opening", progress) as current:
        body_parts = _by_area(opening(device, params.body_min_width / 2))
        current.count(parts=len(body_parts))

    layers = (body_parts, box_parts, handle_parts)
    with stage("bodies", progress, total=sum(map(len, layers))) as current:
        components = layered_components(layers, seed=(2, anchor), advance=current.advance)
        current.count(bodies=len(components))
    grouped: list[list[list[Polygon]]] = [[[], [], []] for _ in components]
    owner = [[0] * len(layer) for layer in layers]
    for index, component in enumerate(components):
        for layer, i in component:
            grouped[index][layer].append(layers[layer][i])
            owner[layer][i] = index

    with stage("candidates", progress) as current:
        candidates = beam_candidates(device, body_parts, params.heal_distance)
        current.count(candidates=len(candidates.geoms))
    with stage("beams", progress, total=len(candidates.geoms)) as current:
        qualified = find_beams(
            candidates,
            body_parts,
            owner[0],
            box_parts,
            margin=margin,
            width_min=params.beam_width_min,
            width_max=width_max,
            thickness=stack.device.thickness,
            net_of=net_of,
            advance=current.advance,
        )
        beams, rejected = qualified.beams, qualified.rejected
        current.count(beams=len(beams), rejected=len(rejected))
    stubs: list[list[Polygon]] = [[] for _ in components]
    for candidate in rejected:
        if len(candidate.bodies) == 1:
            stubs[candidate.bodies[0]].append(candidate.polygon)

    bodies = []
    with stage("mass", progress, total=len(components)) as current:
        for index, group in enumerate(grouped):
            polygons = [MultiPolygon(layer) for layer in group]
            massive = MultiPolygon([*group[0], *stubs[index]])
            inside = openings_within(opening_tree, opening_parts, massive)
            bodies.append(
                Body(
                    index=index,
                    anchor=index == 0,
                    device=polygons[0],
                    box=polygons[1],
                    handle=polygons[2],
                    stubs=MultiPolygon(stubs[index]),
                    inertial=mass_properties(
                        [massive, polygons[1], polygons[2]], stack, openings=inside
                    ),
                    nets=tuple(net_of(part) for part in group[0]),
                )
            )
            current.advance()

    connectors = len(beams) + len(rejected) + len(nets.vias) + len(nets.wires)
    with stage("electrodes", progress, total=connectors) as current:
        electrodes, links = find_electrodes(
            layers,
            owner,
            nets,
            qualified,
            stack.device.material.resistivity,
            params.via_resistance,
            current.advance,
        )
        current.count(electrodes=len(electrodes), links=len(links))
    return Structure(
        bodies=tuple(bodies),
        beams=tuple(beams),
        rejected=tuple(rejected),
        stack=stack,
        params=params,
        nets=nets,
        electrodes=tuple(electrodes),
        links=tuple(links),
        openings=openings,
    )


def _by_area(geometry: MultiPolygon) -> list[Polygon]:
    return sorted(geometry.geoms, key=lambda polygon: polygon.area, reverse=True)

"""Layout to structure."""

from __future__ import annotations

import shapely
from shapely import MultiPolygon, Polygon, STRtree
from shapely.geometry.base import BaseGeometry

from connflex_core.electrical.nets import find_nets
from connflex_core.geometry.beams import find_beams
from connflex_core.geometry.bodies import layered_components, select_anchor
from connflex_core.geometry.heal import as_multipolygon, heal, opening, snap
from connflex_core.geometry.release import fill_release_holes, surviving_box
from connflex_core.layout import Layout
from connflex_core.mechanics.inertia import mass_properties
from connflex_core.model import Body, Structure
from connflex_core.params import ExtractionParams
from connflex_core.stack import Stack


def extract(layout: Layout, stack: Stack, params: ExtractionParams) -> Structure:
    """Extract rigid bodies, beams, and electrical nets from a layout.

    Layers are healed and simplified. Release holes are filled and, with the release layer,
    give the surviving buried oxide. An anchor handle polygon is selected and every polygon
    with no connection to it through the stack is dropped. Bodies are assembled by overlap
    through the opened device layer, oxide, and handle; each gets mass properties from the
    stack with release openings removed from its device layer. Device regions removed by
    the opening are qualified as beams. Remaining device and handle polygons are conductors
    of the nets, joined through vias.
    """
    tolerance = params.simplify_tolerance
    margin = params.beam_overlap_margin
    width_max = params.beam_width_max
    if tolerance is None or margin is None or width_max is None:
        raise ValueError("extraction parameters are not resolved")

    device = heal(layout.device, params.heal_distance, tolerance)
    handle = heal(layout.handle, params.heal_distance, tolerance)
    release = None
    if layout.release is not None:
        release = as_multipolygon(layout.release.simplify(tolerance, preserve_topology=True))
    if params.precision_grid is not None:
        device = snap(device, params.precision_grid)
        handle = snap(handle, params.precision_grid)
        release = None if release is None else snap(release, params.precision_grid)

    device, holes = fill_release_holes(device, params.release_hole_max_area)
    box = surviving_box(
        device, handle, release, holes, params.release_distance, params.buffer_resolution
    )
    openings = shapely.union_all([*holes, *(release.geoms if release is not None else [])])

    handle_parts = _by_area(handle)
    anchor = select_anchor(handle_parts, params.anchor)

    device_parts, box_parts = list(device.geoms), list(box.geoms)
    reachable = layered_components([device_parts, box_parts, handle_parts], seed=(2, anchor))[0]
    kept = [sorted(i for layer, i in reachable if layer == which) for which in range(3)]
    device_parts = [device_parts[i] for i in kept[0]]
    box_parts = [box_parts[i] for i in kept[1]]
    handle_parts = [handle_parts[i] for i in kept[2]]
    anchor = kept[2].index(anchor)
    device = MultiPolygon(device_parts)

    vias = [] if layout.vias is None else list(layout.vias.geoms)
    nets = find_nets(device_parts, handle_parts, vias)
    conductors = STRtree(device_parts)

    def net_of(polygon: Polygon) -> int:
        (i,) = conductors.query(polygon.representative_point(), predicate="within")
        return nets.conductors[int(i)].net

    body_parts = _by_area(opening(device, params.body_min_width / 2))
    components = layered_components([body_parts, box_parts, handle_parts], seed=(2, anchor))
    layers = (body_parts, box_parts, handle_parts)
    bodies = []
    body_of = [0] * len(body_parts)
    for index, component in enumerate(components):
        grouped: list[list[Polygon]] = [[], [], []]
        for layer, i in component:
            grouped[layer].append(layers[layer][i])
            if layer == 0:
                body_of[i] = index
        polygons = [MultiPolygon(group) for group in grouped]
        massive = [_without(polygons[0], openings), polygons[1], polygons[2]]
        bodies.append(
            Body(
                index=index,
                anchor=index == 0,
                device=polygons[0],
                box=polygons[1],
                handle=polygons[2],
                inertial=mass_properties(massive, stack),
                nets=tuple(net_of(part) for part in grouped[0]),
            )
        )

    beams, rejected = find_beams(
        device,
        body_parts,
        body_of,
        margin=margin,
        width_min=params.beam_width_min,
        width_max=width_max,
        thickness=stack.device.thickness,
        net_of=net_of,
    )
    return Structure(
        bodies=tuple(bodies),
        beams=tuple(beams),
        rejected=tuple(rejected),
        material=stack.device.material,
        nets=nets,
    )


def _without(geometry: MultiPolygon, openings: BaseGeometry) -> MultiPolygon:
    return geometry if openings.is_empty else as_multipolygon(geometry.difference(openings))


def _by_area(geometry: MultiPolygon) -> list[Polygon]:
    return sorted(geometry.geoms, key=lambda polygon: polygon.area, reverse=True)

"""Electrodes and links between them."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

from shapely import MultiPolygon, Polygon, STRtree

from connflex_core.electrical.nets import Conductor, Nets

if TYPE_CHECKING:
    from connflex_core.geometry.beams import Qualified


@dataclass(frozen=True)
class Electrode:
    """Silicon of one conductor that moves with one body: a node of the electrical model.

    Attributes:
        index: Position in ``Structure.electrodes``.
        layer: ``"device"`` or ``"handle"``.
        body: Index of the body it moves with.
        conductor: Index of its conductor in ``Nets.conductors``.
        net: Index of its net.
        region: Its silicon, in meters. For a device electrode, one device polygon of the
            body followed by the stubs that touch that polygon alone, or a rejected
            candidate that touches no body polygon, which is an electrode by itself.
    """

    index: int
    layer: str
    body: int
    conductor: int
    net: int
    region: MultiPolygon


@dataclass(frozen=True)
class Link:
    """Conduction between two electrodes.

    Attributes:
        index: Position in ``Structure.links``.
        electrodes: Indices of the two electrodes joined, ascending.
        kind: ``"beam"`` for a beam, ``"bridge"`` for a rejected candidate that touches
            several body polygons, ``"via"`` for a via, ``"wire"`` for a connection
            made off the die.
        source: Index of that beam in ``Structure.beams``, candidate in
            ``Structure.rejected``, via in ``Nets.vias``, or wire in ``Nets.wires``.
        resistance: Resistance in ohms: equation (41) for a beam, ``via_resistance`` for a
            via, zero for a bridge, the resistance of the wire for a wire. None for a beam
            when the device material has no resistivity.
    """

    index: int
    electrodes: tuple[int, int]
    kind: str
    source: int
    resistance: float | None


def find_electrodes(
    layers: Sequence[Sequence[Polygon]],
    owner: Sequence[Sequence[int]],
    nets: Nets,
    qualified: Qualified,
    resistivity: float | None,
    via_resistance: float,
    advance: Callable[[], None] | None = None,
) -> tuple[list[Electrode], list[Link]]:
    """Electrodes of every body and links that conduct between them.

    ``layers`` are device polygons of bodies, oxide polygons, and handle polygons, and
    ``owner`` gives the body of each. Device electrodes come first, then handle
    electrodes, each ordered by body and then by descending area. ``resistivity`` of the
    device layer sets resistance of beam links and ``via_resistance`` that of each via
    link. Each wire of ``nets`` links the electrodes that cover its two ends; an end on
    a beam or a bridge, which belong to no electrode, or a wire with both ends on one
    electrode is an error. ``advance`` is called once per beam, rejected candidate, via,
    and wire.
    """
    parts, oxide, handle = layers
    device = [conductor for conductor in nets.conductors if conductor.layer == "device"]
    tree = STRtree([conductor.polygon for conductor in device])

    def conductor_of(polygon: Polygon) -> Conductor:
        return device[int(tree.nearest(polygon.representative_point()))]

    regions: list[tuple[int, Conductor, list[Polygon]]] = [
        (owner[0][j], conductor_of(part), [part]) for j, part in enumerate(parts)
    ]
    bridges: list[tuple[int, tuple[int, ...]]] = []
    carriers = STRtree(list(oxide))
    for index, touched in enumerate(qualified.rejected_parts):
        polygon = qualified.rejected[index].polygon
        if len(touched) == 1:
            regions[touched[0]][2].append(polygon)
        elif touched:
            bridges.append((index, touched))
        else:
            body = _carrier(polygon, carriers, oxide, owner[1])
            regions.append((body, conductor_of(polygon), [polygon]))
        if advance is not None:
            advance()

    order = sorted(range(len(regions)), key=lambda k: (regions[k][0], -regions[k][2][0].area))
    position = {k: index for index, k in enumerate(order)}
    electrodes = [
        Electrode(index, "device", body, conductor.index, conductor.net, MultiPolygon(polygons))
        for index, (body, conductor, polygons) in enumerate(regions[k] for k in order)
    ]
    on_handle: dict[int, int] = {}
    for i in sorted(range(len(handle)), key=lambda i: (owner[2][i], -handle[i].area)):
        conductor = nets.conductors[len(device) + i]
        on_handle[conductor.index] = len(electrodes)
        electrodes.append(
            Electrode(
                len(electrodes),
                "handle",
                owner[2][i],
                conductor.index,
                conductor.net,
                MultiPolygon([handle[i]]),
            )
        )

    links: list[Link] = []

    def join(a: int, b: int, kind: str, source: int, resistance: float | None = 0.0) -> None:
        links.append(Link(len(links), (min(a, b), max(a, b)), kind, source, resistance))

    for beam, touched in zip(qualified.beams, qualified.beam_parts, strict=True):
        ends = [min(j for j in touched if owner[0][j] == body) for body in beam.bodies]
        resistance = None
        if resistivity is not None:
            resistance = resistivity * beam.length / (beam.width * beam.thickness)
        join(position[ends[0]], position[ends[1]], "beam", beam.index, resistance)
        if advance is not None:
            advance()
    for index, touched in bridges:
        for other in touched[1:]:
            join(position[touched[0]], position[other], "bridge", index)
    for via in nets.vias:
        if advance is not None:
            advance()
        for reached in via.device:
            above = max(
                (e for e in electrodes if e.conductor == reached),
                key=lambda e: (
                    via.polygon.intersection(e.region).area,
                    -via.polygon.distance(e.region),
                ),
                default=None,
            )
            if above is not None:
                for below in via.handle:
                    join(above.index, on_handle[below], "via", via.index, via_resistance)
    for connection in nets.wires:
        ends = []
        for (point, layer), under in zip(connection.wire.ends, connection.conductors, strict=True):
            if layer == "handle":
                ends.append(on_handle[under])
                continue
            covering = [
                e.index for e in electrodes if e.conductor == under and e.region.covers(point)
            ]
            if not covering:
                raise ValueError(
                    f"wire {connection.index}: end at ({point.x:g}, {point.y:g}) m lies on a beam "
                    f"or bridge, which belongs to no electrode; place it on a pad"
                )
            ends.append(covering[0])
        if ends[0] == ends[1]:
            raise ValueError(f"wire {connection.index}: both ends lie on one electrode, {ends[0]}")
        join(ends[0], ends[1], "wire", connection.index, connection.wire.resistance)
        if advance is not None:
            advance()
    return electrodes, links


def _carrier(
    polygon: Polygon, tree: STRtree, oxide: Sequence[Polygon], owner: Sequence[int]
) -> int:
    hits = tree.query(polygon, predicate="intersects")
    if hits.size == 0:
        return 0
    return owner[int(max(hits, key=lambda i: polygon.intersection(oxide[i]).area))]

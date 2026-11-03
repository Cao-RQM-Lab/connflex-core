"""Extracted structure."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import networkx as nx
from shapely import MultiPolygon, Polygon
from shapely.geometry import mapping

from connflex_core.mechanics.inertia import MassProperties


@dataclass(frozen=True)
class Body:
    """A rigid body and the polygons it occupies in each layer, in meters.

    Attributes:
        index: Position in ``Structure.bodies``. The anchor is 0.
        anchor: Whether this is the immovable anchor body.
        device: Device-layer polygons. Empty when the body has none in this layer.
        box: Surviving buried oxide polygons.
        handle: Handle-layer polygons.
        inertial: Mass, center of mass, and inertia tensor.
    """

    index: int
    anchor: bool
    device: MultiPolygon
    box: MultiPolygon
    handle: MultiPolygon
    inertial: MassProperties


@dataclass(frozen=True)
class Beam:
    """A straight beam joining two bodies.

    Attributes:
        index: Position in ``Structure.beams``.
        bodies: Indices of the two bodies joined, ascending.
        length: Distance between the attachment points, in meters.
        width: In-plane width: the beam area divided by the length.
        thickness: Device layer thickness.
        start: Attachment point on the first body, (x, y, z) at the device mid-plane.
        end: Attachment point on the second body.
        polygon: The device-layer region of the beam.
    """

    index: int
    bodies: tuple[int, int]
    length: float
    width: float
    thickness: float
    start: tuple[float, float, float]
    end: tuple[float, float, float]
    polygon: Polygon


@dataclass(frozen=True)
class RejectedCandidate:
    """A device-layer region narrower than a body that did not qualify as a beam.

    Attributes:
        polygon: The region.
        reason: Why it was rejected.
    """

    polygon: Polygon
    reason: str


@dataclass(frozen=True)
class Structure:
    """Rigid bodies and beams extracted from a layout.

    Attributes:
        bodies: All bodies, the anchor first.
        beams: Accepted beams.
        rejected: Beam candidates that did not qualify, with the reason for each.
    """

    bodies: tuple[Body, ...]
    beams: tuple[Beam, ...]
    rejected: tuple[RejectedCandidate, ...]

    @property
    def anchor(self) -> Body:
        """The immovable anchor body."""
        return self.bodies[0]

    @property
    def graph(self) -> nx.MultiGraph[int]:
        """Bodies as nodes and beams as edges, with the objects as ``body`` and ``beam`` attributes."""
        graph: nx.MultiGraph[int] = nx.MultiGraph()
        for body in self.bodies:
            graph.add_node(body.index, body=body)
        for beam in self.beams:
            graph.add_edge(*beam.bodies, key=beam.index, beam=beam)
        return graph

    def to_dict(self) -> dict[str, Any]:
        """JSON-serializable tree with polygons as GeoJSON-style mappings."""
        return {
            "bodies": [
                {
                    "index": body.index,
                    "anchor": body.anchor,
                    "device": mapping(body.device),
                    "box": mapping(body.box),
                    "handle": mapping(body.handle),
                    "mass": body.inertial.mass,
                    "center": body.inertial.center.tolist(),
                    "inertia": body.inertial.inertia.tolist(),
                }
                for body in self.bodies
            ],
            "beams": [
                {
                    "index": beam.index,
                    "bodies": list(beam.bodies),
                    "length": beam.length,
                    "width": beam.width,
                    "thickness": beam.thickness,
                    "start": list(beam.start),
                    "end": list(beam.end),
                    "polygon": mapping(beam.polygon),
                }
                for beam in self.beams
            ],
            "rejected": [
                {"polygon": mapping(candidate.polygon), "reason": candidate.reason}
                for candidate in self.rejected
            ],
        }

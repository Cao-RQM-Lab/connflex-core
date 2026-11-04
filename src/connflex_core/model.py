"""Extracted structure."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import networkx as nx
import numpy as np
from shapely import MultiPolygon, Polygon
from shapely.geometry import mapping

from connflex_core.electrical.nets import Nets
from connflex_core.materials import Material
from connflex_core.mechanics.inertia import MassProperties
from connflex_core.params import SolveOptions
from connflex_core.progress import Progress

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from connflex_core.mechanics.modal import ModalResult


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
        nets: Net of each polygon in ``device``, in order.
    """

    index: int
    anchor: bool
    device: MultiPolygon
    box: MultiPolygon
    handle: MultiPolygon
    inertial: MassProperties
    nets: tuple[int, ...]


@dataclass(frozen=True)
class Beam:
    """A straight beam joining two bodies.

    Attributes:
        index: Position in ``Structure.beams``.
        bodies: Indices of the two bodies joined, ascending.
        length: Distance between attachment points, in meters.
        width: In-plane width: beam area divided by length.
        thickness: Device layer thickness.
        start: Attachment point on the first body, (x, y, z) at device mid-plane.
        end: Attachment point on the second body.
        polygon: Device-layer region of the beam.
        net: Net of the device conductor the beam belongs to.
    """

    index: int
    bodies: tuple[int, int]
    length: float
    width: float
    thickness: float
    start: tuple[float, float, float]
    end: tuple[float, float, float]
    polygon: Polygon
    net: int


@dataclass(frozen=True)
class RejectedCandidate:
    """A device-layer region narrower than a body that did not qualify as a beam.

    Attributes:
        polygon: Rejected region.
        reason: Why it was rejected.
        bodies: Indices of the bodies it touches, ascending.
    """

    polygon: Polygon
    reason: str
    bodies: tuple[int, ...]


@dataclass(frozen=True)
class Structure:
    """Rigid bodies and beams extracted from a layout.

    Attributes:
        bodies: All bodies, anchor first.
        beams: Accepted beams.
        rejected: Beam candidates that did not qualify, each with a reason.
        material: Device-layer material, which sets beam stiffness.
        nets: Electrical nets of conductors that remain after pruning.
    """

    bodies: tuple[Body, ...]
    beams: tuple[Beam, ...]
    rejected: tuple[RejectedCandidate, ...]
    material: Material
    nets: Nets

    @property
    def anchor(self) -> Body:
        """Immovable anchor body."""
        return self.bodies[0]

    @property
    def graph(self) -> nx.MultiGraph[int]:
        """Bodies as nodes and beams as edges, with objects as ``body`` and ``beam`` attributes."""
        graph: nx.MultiGraph[int] = nx.MultiGraph()
        for body in self.bodies:
            graph.add_node(body.index, body=body)
        for beam in self.beams:
            graph.add_edge(*beam.bodies, key=beam.index, beam=beam)
        return graph

    def modal(
        self, options: SolveOptions | None = None, *, progress: Progress | None = None
    ) -> ModalResult:
        """Resonant frequencies and mode shapes of the assembled system."""
        from connflex_core.mechanics.modal import modal

        return modal(self, options, progress=progress)

    def static(
        self, loads: Mapping[int, Sequence[float]], *, progress: Progress | None = None
    ) -> NDArray[np.float64]:
        """Displacements under static generalized forces, keyed by body index."""
        from connflex_core.mechanics.modal import static

        return static(self, loads, progress=progress)

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
                    "nets": list(body.nets),
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
                    "net": beam.net,
                }
                for beam in self.beams
            ],
            "rejected": [
                {
                    "polygon": mapping(candidate.polygon),
                    "reason": candidate.reason,
                    "bodies": list(candidate.bodies),
                }
                for candidate in self.rejected
            ],
            "nets": self.nets.to_dict(),
        }

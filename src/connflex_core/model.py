"""Extracted structure."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import networkx as nx
import numpy as np
from shapely import MultiPolygon, Point, Polygon
from shapely.geometry import mapping

from connflex_core.electrical.electrodes import Electrode, Link
from connflex_core.electrical.nets import Nets
from connflex_core.mechanics.inertia import MassProperties
from connflex_core.params import ElectrostaticParams, ExtractionParams, SolveOptions, StaticOptions
from connflex_core.progress import Progress
from connflex_core.stack import Stack

if TYPE_CHECKING:
    from numpy.typing import ArrayLike, NDArray

    from connflex_core.coupled.static import Equilibrium
    from connflex_core.electrostatics.capacitance import CapacitanceModel
    from connflex_core.electrostatics.plates import Plate
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
        stubs: Device-layer regions attached to this body alone: beam candidates that touch
            no other body. They move with the body and count in ``inertial``, but they are
            not part of ``device``, which is the body opening that defines connectivity and
            beam qualification.
        inertial: Mass, center of mass, and inertia tensor over device, stubs, box, and
            handle, with release openings removed.
        nets: Net of each polygon in ``device``, in order.
    """

    index: int
    anchor: bool
    device: MultiPolygon
    box: MultiPolygon
    handle: MultiPolygon
    stubs: MultiPolygon
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
        stack: Layer stack of the extraction. Its device material sets beam stiffness.
        params: Extraction parameters it was extracted with.
        nets: Electrical nets of conductors that remain after pruning.
        electrodes: Nodes of the electrical model, device electrodes first, then handle
            electrodes, each ordered by body and then by descending area.
        links: Beams, bridges, vias, and wires that conduct between electrodes.
        openings: Release openings as disjoint polygons: release holes that were filled in
            device polygons, and polygons of the release layer. Device silicon is absent
            there, although ``Body.device`` and ``Electrode.region`` cover them.
    """

    bodies: tuple[Body, ...]
    beams: tuple[Beam, ...]
    rejected: tuple[RejectedCandidate, ...]
    stack: Stack
    params: ExtractionParams
    nets: Nets
    electrodes: tuple[Electrode, ...]
    links: tuple[Link, ...]
    openings: MultiPolygon

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

    def electrode_at(self, point: Point | tuple[float, float], layer: str) -> Electrode:
        """Electrode of ``layer`` covering ``point``."""
        point = Point(point)
        for electrode in self.electrodes:
            if electrode.layer == layer and electrode.region.covers(point):
                return electrode
        raise ValueError(f"no {layer} electrode at ({point.x}, {point.y})")

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

    def plates(
        self, params: ElectrostaticParams | None = None, *, progress: Progress | None = None
    ) -> tuple[Plate, ...]:
        """Capacitors between device electrodes and handle electrodes beneath them."""
        from connflex_core.electrostatics.plates import find_plates

        return find_plates(self, params, progress=progress)

    def capacitance(
        self, params: ElectrostaticParams | None = None, *, progress: Progress | None = None
    ) -> CapacitanceModel:
        """Capacitance between all electrodes at rest, from plates and sidewalls."""
        from connflex_core.electrostatics.capacitance import find_capacitance

        return find_capacitance(self, params, progress=progress)

    def equilibrium(
        self,
        model: CapacitanceModel,
        potentials: Mapping[int, float],
        *,
        charges: Mapping[int, float] | None = None,
        loads: Mapping[int, Sequence[float]] | None = None,
        start: Equilibrium | ArrayLike | None = None,
        options: StaticOptions | None = None,
        progress: Progress | None = None,
    ) -> Equilibrium:
        """Equilibrium under potentials and charges on nets, found from rest or from
        ``start``, with the capacitance model of this structure."""
        from connflex_core.coupled.static import equilibrium

        return equilibrium(
            self,
            model,
            potentials,
            charges=charges,
            loads=loads,
            start=start,
            options=options,
            progress=progress,
        )

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
                    "stubs": mapping(body.stubs),
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
            "electrodes": [
                {
                    "index": electrode.index,
                    "layer": electrode.layer,
                    "body": electrode.body,
                    "conductor": electrode.conductor,
                    "net": electrode.net,
                    "region": mapping(electrode.region),
                }
                for electrode in self.electrodes
            ],
            "links": [
                {
                    "index": link.index,
                    "electrodes": list(link.electrodes),
                    "kind": link.kind,
                    "source": link.source,
                    "resistance": link.resistance,
                }
                for link in self.links
            ],
        }

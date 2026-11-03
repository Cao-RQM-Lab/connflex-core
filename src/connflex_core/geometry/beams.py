"""Beam qualification."""

from __future__ import annotations

from collections.abc import Callable, Sequence

import shapely
from shapely import MultiPolygon, Polygon, STRtree

from connflex_core.geometry.heal import as_multipolygon
from connflex_core.model import Beam, RejectedCandidate


def find_beams(
    device: MultiPolygon,
    body_parts: Sequence[Polygon],
    body_of: Sequence[int],
    *,
    margin: float,
    width_min: float,
    width_max: float,
    thickness: float,
    net_of: Callable[[Polygon], int],
) -> tuple[list[Beam], list[RejectedCandidate]]:
    """Qualify device regions outside the bodies as beams joining exactly two bodies.

    ``body_parts`` are the device-layer polygons of all bodies and ``body_of`` maps each part
    to its body index. ``net_of`` gives the net of a device-layer region.
    """
    tree = STRtree(list(body_parts))
    candidates = as_multipolygon(device.difference(shapely.union_all(list(body_parts))))
    z = -thickness / 2
    beams: list[Beam] = []
    rejected: list[RejectedCandidate] = []
    for candidate in candidates.geoms:
        expanded = candidate.buffer(margin, join_style="mitre")
        touched: dict[int, list[Polygon]] = {}
        for j in tree.query(expanded, predicate="intersects"):
            if expanded.intersection(body_parts[j]).area > 0:
                touched.setdefault(body_of[j], []).append(body_parts[j])
        if len(touched) != 2:
            noun = "body" if len(touched) == 1 else "bodies"
            rejected.append(RejectedCandidate(candidate, f"touches {len(touched)} {noun}"))
            continue
        (a, parts_a), (b, parts_b) = sorted(touched.items())
        start = expanded.intersection(shapely.union_all(parts_a)).centroid
        end = expanded.intersection(shapely.union_all(parts_b)).centroid
        length = start.distance(end)
        if length <= 0:
            rejected.append(RejectedCandidate(candidate, "zero length"))
            continue
        width = candidate.area / length
        if width < width_min:
            reason = f"width {width:.3g} m is below beam_width_min {width_min:.3g} m"
            rejected.append(RejectedCandidate(candidate, reason))
            continue
        if width > width_max:
            reason = f"width {width:.3g} m is above beam_width_max {width_max:.3g} m"
            rejected.append(RejectedCandidate(candidate, reason))
            continue
        beams.append(
            Beam(
                index=len(beams),
                bodies=(a, b),
                length=length,
                width=width,
                thickness=thickness,
                start=(start.x, start.y, z),
                end=(end.x, end.y, z),
                polygon=candidate,
                net=net_of(candidate),
            )
        )
    return beams, rejected

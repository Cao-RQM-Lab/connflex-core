"""Beam qualification."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import NamedTuple

import shapely
from shapely import MultiPolygon, Polygon, STRtree

from connflex_core.geometry.heal import as_multipolygon, opening
from connflex_core.model import Beam, RejectedCandidate


def beam_candidates(
    device: MultiPolygon, body_parts: Sequence[Polygon], sliver: float
) -> MultiPolygon:
    """Device regions outside the bodies, with parts narrower than ``sliver`` removed.

    The body opening does not reproduce device edges bit for bit, so the difference leaves
    hairlines along them that would join separate candidates into one.
    """
    outside = as_multipolygon(device.difference(shapely.union_all(list(body_parts))))
    return opening(outside, sliver / 2)


class Qualified(NamedTuple):
    """Beams and rejected candidates, each with indices of the body polygons it touches."""

    beams: list[Beam]
    rejected: list[RejectedCandidate]
    beam_parts: list[tuple[int, ...]]
    rejected_parts: list[tuple[int, ...]]


def find_beams(
    candidates: MultiPolygon,
    body_parts: Sequence[Polygon],
    body_of: Sequence[int],
    *,
    margin: float,
    width_min: float,
    width_max: float,
    thickness: float,
    net_of: Callable[[Polygon], int],
    advance: Callable[[], None] | None = None,
) -> Qualified:
    """Qualify beam candidates as beams joining exactly two bodies.

    ``body_parts`` are the device-layer polygons of all bodies and ``body_of`` maps each part
    to its body index. ``net_of`` gives the net of a device-layer region. ``advance`` is
    called once per candidate.
    """
    tree = STRtree(list(body_parts))
    z = -thickness / 2
    found = Qualified([], [], [], [])
    beams = found.beams

    def reject(candidate: Polygon, reason: str, bodies: tuple[int, ...]) -> None:
        found.rejected.append(RejectedCandidate(candidate, reason, bodies))
        found.rejected_parts.append(parts)

    for candidate in candidates.geoms:
        if advance is not None:
            advance()
        expanded = candidate.buffer(margin, join_style="mitre")
        touched: dict[int, list[int]] = {}
        for j in tree.query(expanded, predicate="intersects"):
            if expanded.intersection(body_parts[j]).area > 0:
                touched.setdefault(body_of[j], []).append(int(j))
        parts = tuple(sorted(j for group in touched.values() for j in group))
        if len(touched) != 2:
            noun = "body" if len(touched) == 1 else "bodies"
            reject(candidate, f"touches {len(touched)} {noun}", tuple(sorted(touched)))
            continue
        (a, parts_a), (b, parts_b) = sorted(touched.items())
        start = expanded.intersection(shapely.union_all([body_parts[j] for j in parts_a])).centroid
        end = expanded.intersection(shapely.union_all([body_parts[j] for j in parts_b])).centroid
        length = start.distance(end)
        if length <= 0:
            reject(candidate, "zero length", (a, b))
            continue
        width = candidate.area / length
        if width < width_min:
            reason = f"width {width:.3g} m is below beam_width_min {width_min:.3g} m"
            reject(candidate, reason, (a, b))
            continue
        if width > width_max:
            reason = f"width {width:.3g} m is above beam_width_max {width_max:.3g} m"
            reject(candidate, reason, (a, b))
            continue
        found.beam_parts.append(parts)
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
    return found

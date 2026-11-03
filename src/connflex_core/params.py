"""Extraction and solve parameters."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from connflex_core._validation import require_nonnegative, require_positive

SOLVE_METHODS = ("dense", "sparse")


@dataclass(frozen=True)
class ExtractionParams:
    """Geometry extraction parameters. Lengths in meters, areas in square meters.

    Attributes:
        release_hole_max_area: Holes inside device polygons at or below this area are release
            holes. They enter the release estimate and are then filled.
        release_distance: Isotropic undercut applied to the surviving buried oxide.
        body_min_width: Features narrower than this are not rigid bodies and become beam
            candidates. Also the default upper width bound of a beam.
        heal_distance: Buffer out then in on the device and handle layers to close hairline gaps.
        simplify_tolerance: Douglas-Peucker tolerance applied after healing. None means
            heal_distance.
        beam_overlap_margin: Expansion of beam candidates so they overlap the bodies they touch.
            None means heal_distance.
        beam_width_min: Candidates with computed width below this are discarded as slivers.
        beam_width_max: Candidates with computed width above this are discarded as non-slender.
            None means body_min_width.
        anchor: None selects the largest handle polygon, an int selects a handle polygon by
            index, and a point (x, y) in meters selects the handle polygon containing it.
        buffer_resolution: Segments per quarter circle for round buffers.
        precision_grid: Snap grid applied before overlay operations. None disables snapping.
    """

    release_hole_max_area: float
    release_distance: float
    body_min_width: float
    heal_distance: float = 10e-9
    simplify_tolerance: float | None = None
    beam_overlap_margin: float | None = None
    beam_width_min: float = 0.0
    beam_width_max: float | None = None
    anchor: int | tuple[float, float] | None = None
    buffer_resolution: int = 8
    precision_grid: float | None = None

    def __post_init__(self) -> None:
        require_positive("release_hole_max_area", self.release_hole_max_area)
        require_positive("release_distance", self.release_distance)
        require_positive("body_min_width", self.body_min_width)
        require_positive("heal_distance", self.heal_distance)
        if self.simplify_tolerance is None:
            object.__setattr__(self, "simplify_tolerance", self.heal_distance)
        else:
            require_positive("simplify_tolerance", self.simplify_tolerance)
        margin = (
            self.heal_distance if self.beam_overlap_margin is None else self.beam_overlap_margin
        )
        require_positive("beam_overlap_margin", margin)
        if not margin < self.body_min_width / 2:
            raise ValueError(
                f"beam_overlap_margin must be smaller than body_min_width / 2, "
                f"got {margin} with body_min_width {self.body_min_width}"
            )
        object.__setattr__(self, "beam_overlap_margin", margin)
        require_nonnegative("beam_width_min", self.beam_width_min)
        width_max = self.body_min_width if self.beam_width_max is None else self.beam_width_max
        require_positive("beam_width_max", width_max)
        if not self.beam_width_min < width_max:
            raise ValueError(
                f"beam_width_min must be smaller than beam_width_max, "
                f"got {self.beam_width_min} and {width_max}"
            )
        object.__setattr__(self, "beam_width_max", width_max)
        if self.buffer_resolution < 1:
            raise ValueError(f"buffer_resolution must be at least 1, got {self.buffer_resolution}")
        if self.precision_grid is not None:
            require_positive("precision_grid", self.precision_grid)
        anchor = self.anchor
        if isinstance(anchor, int):
            require_nonnegative("anchor", anchor)
        elif anchor is not None and (len(anchor) != 2):
            raise ValueError(
                f"anchor must be None, a handle polygon index, or an (x, y) point, got {anchor!r}"
            )


@dataclass(frozen=True)
class SolveOptions:
    """Modal solve options.

    Attributes:
        n_modes: Number of lowest modes to return. None returns all modes.
        method: "dense" solves the full generalized eigenproblem; "sparse" uses shift-invert
            iteration.
        sigma: Shift for the sparse solver in squared angular frequency. Ignored by the dense
            solver.
    """

    n_modes: int | None = None
    method: Literal["dense", "sparse"] = "dense"
    sigma: float = 0.0

    def __post_init__(self) -> None:
        if self.n_modes is not None and self.n_modes < 1:
            raise ValueError(f"n_modes must be at least 1, got {self.n_modes}")
        if self.method not in SOLVE_METHODS:
            raise ValueError(f"method must be one of {SOLVE_METHODS}, got {self.method!r}")

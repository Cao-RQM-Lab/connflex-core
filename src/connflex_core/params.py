"""Extraction and solve parameters."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from connflex_core._validation import require_at_least, require_nonnegative, require_positive

SOLVE_METHODS = ("dense", "sparse")
SIDEWALL_METHODS = ("field", "plates")


@dataclass(frozen=True)
class ExtractionParams:
    """Extraction parameters. Lengths in meters, areas in square meters, resistances in ohms.

    Attributes:
        release_hole_max_area: Holes inside device polygons at or below this area are release
            holes. They enter the release estimate, are removed from device-layer mass, and
            are filled for connectivity and beam extraction.
        release_distance: Isotropic undercut applied to the surviving buried oxide.
        body_min_width: Features narrower than this are not rigid bodies and become beam
            candidates. Also the default upper width bound of a beam.
        heal_distance: Buffer out then in on device and handle layers to close hairline gaps.
            Beam candidates narrower than this are dropped as slivers of the body opening.
        simplify_tolerance: Douglas-Peucker tolerance applied after healing. None means
            heal_distance.
        beam_overlap_margin: Expansion of beam candidates so they overlap the bodies they touch.
            None means heal_distance.
        beam_width_min: Candidates with computed width below this are discarded as slivers.
        beam_width_max: Candidates with computed width above this are discarded as non-slender.
            None means body_min_width.
        anchor: None selects the largest handle polygon, a point (x, y) in meters selects the
            handle polygon containing it, and an int selects a handle polygon by index in
            order of descending area.
        buffer_resolution: Segments per quarter circle of the round release undercut. Healing
            and body opening use mitered corners.
        precision_grid: Snap grid applied before overlay operations. None disables snapping.
        via_resistance: Resistance of each via link between a device and a handle electrode,
            whatever the via's area. Zero is an ideal connection.
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
    via_resistance: float = 0.0

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
        require_nonnegative("via_resistance", self.via_resistance)
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


@dataclass(frozen=True)
class StaticOptions:
    """Options of the static equilibrium under electrode potentials.

    Attributes:
        tolerance: Newton iteration stops when the strain energy of its last step is
            below this fraction, squared, of the strain energy of the displacement, or
            one step after it is below this fraction itself.
        max_iterations: Most Newton iterations at one set of potentials.
        steps: Equal steps in which potentials rise from zero to the values asked for.
            The equilibrium is followed from one to the next, and each is kept in the
            result. A step that finds no stable equilibrium is halved.
        pull_in_resolution: Halving stops, and pull-in is located, when a step is
            smaller than this fraction of the potentials asked for.
    """

    tolerance: float = 1e-9
    max_iterations: int = 30
    steps: int = 1
    pull_in_resolution: float = 1e-3

    def __post_init__(self) -> None:
        require_positive("tolerance", self.tolerance)
        require_at_least("max_iterations", self.max_iterations, 1)
        require_at_least("steps", self.steps, 1)
        if not 0 < self.pull_in_resolution < 1:
            raise ValueError(
                f"pull_in_resolution must be between 0 and 1, got {self.pull_in_resolution}"
            )


@dataclass(frozen=True)
class ElectrostaticParams:
    """Electrostatic model parameters. Lengths in meters.

    Attributes:
        ambient_permittivity: Relative permittivity of the medium that fills released gaps
            and trenches.
        element_size: Longest boundary element on conductor outlines. Sidewall capacitance
            requires it; None serves analyses that mesh no outline.
        corner_grading: Ratio of ``element_size`` to the element length asked for at a
            convex corner of a right angle or sharper. The length asked for grows with
            distance from the corner. One disables refinement.
        max_elements: Most boundary elements solved as one dense system. A larger mesh is
            solved in windows, none of which grows past this number.
        window_elements: Most boundary elements in the core of a window, the part of a
            large mesh whose charge one dense solution supplies.
        window_tolerance: Overlap of a window around its core doubles until the charge
            on the core changes by less than this fraction.
        sidewall_method: How sidewall capacitance is found: ``"plates"`` takes facing
            walls as parallel plates, ``"field"`` solves the field between conductors.
            No default: which one suits a device is the caller's judgment of its
            geometry, and sidewall capacitance requires the choice.
        interaction_distance: Device electrodes whose silicon lies within this distance
            of each other are solved together for sidewall capacitance, as are chains of
            such neighbors; electrodes farther apart share none. The capacitance model
            of a structure requires it when sidewalls are included. With method
            ``"plates"`` it is also the widest gap across which two walls count, and
            with method ``"field"`` the widest overlap of a window.
        include_plates: Whether the capacitance model of a structure holds plates between
            device and handle electrodes.
        include_sidewalls: Whether it holds sidewall capacitance between device
            electrodes.
    """

    ambient_permittivity: float = 1.0
    element_size: float | None = None
    corner_grading: float = 16.0
    max_elements: int = 8000
    window_elements: int = 500
    window_tolerance: float = 1e-3
    sidewall_method: Literal["field", "plates"] | None = None
    interaction_distance: float | None = None
    include_plates: bool = True
    include_sidewalls: bool = True

    def __post_init__(self) -> None:
        require_at_least("ambient_permittivity", self.ambient_permittivity, 1)
        if self.element_size is not None:
            require_positive("element_size", self.element_size)
        require_at_least("corner_grading", self.corner_grading, 1)
        require_at_least("max_elements", self.max_elements, 1)
        require_at_least("window_elements", self.window_elements, 1)
        require_positive("window_tolerance", self.window_tolerance)
        if self.sidewall_method not in (None, *SIDEWALL_METHODS):
            raise ValueError(
                f"sidewall_method must be one of {SIDEWALL_METHODS}, got {self.sidewall_method!r}"
            )
        if self.interaction_distance is not None:
            require_positive("interaction_distance", self.interaction_distance)

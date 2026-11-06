"""Extraction and solve parameters."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from connflex_core._validation import require_at_least, require_nonnegative, require_positive

SOLVE_METHODS = ("dense", "sparse")


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
        max_iterations: Most Newton iterations. Reaching it raises ``NoEquilibriumError``.
    """

    tolerance: float = 1e-9
    max_iterations: int = 30

    def __post_init__(self) -> None:
        require_positive("tolerance", self.tolerance)
        require_at_least("max_iterations", self.max_iterations, 1)


@dataclass(frozen=True)
class SmallSignalOptions:
    """Options of the small-signal model about an equilibrium.

    Attributes:
        quality_factor: Quality factor of every mode at the equilibrium, equation (132).
            Required: damping is an input, not computed.
        linear_tolerance: ``SmallSignal.linearity`` warns when either of its ratios,
            equation (137), exceeds this.
    """

    quality_factor: float
    linear_tolerance: float = 0.1

    def __post_init__(self) -> None:
        require_positive("quality_factor", self.quality_factor)
        require_positive("linear_tolerance", self.linear_tolerance)


@dataclass(frozen=True)
class TransientOptions:
    """Options of a transient run.

    Attributes:
        quality_factor: Quality factor of every mode of the structure without
            potentials, which sets its damping, equation (140). Required.
        tolerance: Error the integrator allows per step, relative to each state and to
            the scales of equation (146).
        max_step: Longest time step in seconds. None leaves step length to error control.
        contact_clearance: The run ends with a contact event when the clearance of
            equation (126) falls to this.
    """

    quality_factor: float
    tolerance: float = 1e-6
    max_step: float | None = None
    contact_clearance: float = 0.1

    def __post_init__(self) -> None:
        require_positive("quality_factor", self.quality_factor)
        require_positive("tolerance", self.tolerance)
        if self.max_step is not None:
            require_positive("max_step", self.max_step)
        if not 0 < self.contact_clearance < 1:
            raise ValueError(
                f"contact_clearance must lie between 0 and 1, got {self.contact_clearance}"
            )


@dataclass(frozen=True)
class ElectrostaticParams:
    """Electrostatic model parameters. Lengths in meters.

    Attributes:
        ambient_permittivity: Relative permittivity of the medium that fills released gaps
            and trenches.
        element_size: Longest boundary element on conductor outlines. Required when a pair
            of nets is solved by the field method; None serves analyses that mesh no
            outline.
        corner_grading: Ratio of ``element_size`` to the element length asked for at a
            convex corner of a right angle or sharper. The length asked for grows with
            distance from the corner. One disables refinement.
        max_elements: Most boundary elements solved as one dense system. A dense problem
            that needs more is solved in windows, none of which grows past this number.
        window_elements: Most boundary elements in the core of a window, the part of a
            large mesh whose charge one dense solution supplies.
        window_tolerance: Overlap of a window around its core doubles until the charge
            on the core changes by less than this fraction.
        dense_cutoff: Pairs of nets whose feature size, the length of their facing walls
            in units of the gap between them, is below this are solved by the field
            method in one dense system.
        plates_cutoff: Pairs whose feature size is at least this take their facing walls
            as parallel plates; pairs between the two cutoffs are solved by the field
            method in windows. At least ``dense_cutoff``. Zero for both puts every pair
            on plates, infinity for both every pair in a dense system.
        interaction_distance: Conductors whose silicon lies within this distance of each
            other interact, and chains of such neighbors form one cluster; conductors
            farther apart share no sidewall capacitance. It is also the widest gap
            across which two walls count as plates and the widest overlap of a window.
            Sidewall capacitance requires it.
        include_plates: Whether the capacitance model of a structure holds plates between
            device and handle electrodes.
        include_sidewalls: Whether it holds sidewall capacitance between device
            electrodes.
        surrogate_tolerance: With bodies displaced, a cluster solved by the field method
            reuses the nearest earlier solution of it, expanded to second order, while
            the error of that expansion in capacitance and its gradient stays below this
            fraction, equations (148) to (150). None solves every displacement anew.
    """

    ambient_permittivity: float = 1.0
    element_size: float | None = None
    corner_grading: float = 16.0
    max_elements: int = 8000
    window_elements: int = 500
    window_tolerance: float = 1e-3
    dense_cutoff: float = 10.0
    plates_cutoff: float = 100.0
    interaction_distance: float | None = None
    include_plates: bool = True
    include_sidewalls: bool = True
    surrogate_tolerance: float | None = None

    def __post_init__(self) -> None:
        require_at_least("ambient_permittivity", self.ambient_permittivity, 1)
        if self.surrogate_tolerance is not None and not 0 < self.surrogate_tolerance < 1:
            raise ValueError(
                f"surrogate_tolerance must lie between 0 and 1, got {self.surrogate_tolerance}"
            )
        if self.element_size is not None:
            require_positive("element_size", self.element_size)
        require_at_least("corner_grading", self.corner_grading, 1)
        require_at_least("max_elements", self.max_elements, 1)
        require_at_least("window_elements", self.window_elements, 1)
        require_positive("window_tolerance", self.window_tolerance)
        if not self.dense_cutoff >= 0:
            raise ValueError(f"dense_cutoff must be at least 0, got {self.dense_cutoff}")
        if not self.plates_cutoff >= self.dense_cutoff:
            raise ValueError(
                f"plates_cutoff must be at least dense_cutoff = {self.dense_cutoff}, "
                f"got {self.plates_cutoff}"
            )
        if self.interaction_distance is not None:
            require_positive("interaction_distance", self.interaction_distance)

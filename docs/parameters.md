# Parameters

Every tunable is a named field on a parameter object or an input object.
There are no hidden constants. A test checks that this page and the fields agree. Each
object is built in Python from keyword arguments and passed to the function
that uses it; the library reads no settings files.

Lengths are in meters, areas in square meters, densities in kilograms per
cubic meter, moduli in pascals, resistivities in ohm meters, resistances in
ohms. Relative permittivities have no unit. Fields marked
**required** have no default:
a wrong release distance or body width yields a plausible but wrong
structure, so the library does not guess them.

## `ExtractionParams`

Symbols refer to the [theory of operation](theory.md).

| Field | Symbol | Default | Effect |
| --- | --- | --- | --- |
| `heal_distance` | $r_h$ | 10 nm | Buffer out then in on device and handle layers to close hairline gaps, equation (1). Beam candidates narrower than $r_h$ are dropped as slivers of the body opening. |
| `simplify_tolerance` | $r_h$ | = `heal_distance` | Douglas–Peucker tolerance applied after healing. |
| `release_hole_max_area` | $A_\text{min}$ | **required** | Holes inside device polygons at or below this area are release holes: they enter the release estimate, are removed from device-layer mass, and are filled for connectivity and beam extraction. Holes are identified before healing and once more after it. |
| `release_distance` | $d$ | **required** | Isotropic undercut of the surviving buried oxide, equation (2). |
| `body_min_width` | $w_\text{min}$ | **required** | Features narrower than this become beam candidates, equation (3). Also the upper width bound of a beam. |
| `beam_overlap_margin` | $\epsilon$ | = `heal_distance` | Expansion of beam candidates so they overlap the bodies they touch, equation (4). |
| `beam_width_min` | | 0 (disabled) | Candidates narrower than this are discarded as slivers. |
| `beam_width_max` | | = `body_min_width` | Candidates wider than this are discarded as non-slender. |
| `anchor` | | largest handle polygon | Override by a point inside the intended handle polygon, or by index into the handle polygons ordered by descending area. |
| `buffer_resolution` | | 8 | Segments per quarter circle of the round release undercut, equation (2). Healing and body opening use mitered corners. |
| `precision_grid` | | off | Optional snap grid applied before overlay operations. |
| `via_resistance` | $R_v$ | 0 | Resistance of each via link between a device and a handle electrode, whatever the via's area. Zero is an ideal connection. |

## `Stack`, `Layer`, and `Material`

A stack lists the SOI layers from top down: `device`, `box`, `handle`. Beam
stiffness uses the device-layer material.

| Field | Symbol | Default | Effect |
| --- | --- | --- | --- |
| `Layer.thickness` | $t_l$ | **required** | Layer thickness, equations (5) and (8). |
| `Material.density` | $\rho_l$ | silicon 2330, oxide 2200 | Mass density, equation (5). |
| `Material.youngs_modulus` | $E$ | 169 GPa for silicon | Young's modulus. Silicon's default is the in-plane $\langle 110 \rangle$ value; beams along $\langle 100 \rangle$ see about 130 GPa. |
| `Material.poisson_ratio` | $\nu$ | 0.28 | Used to derive $G = E / (2(1 + \nu))$ when `shear_modulus` is not given. |
| `Material.shear_modulus` | $G$ | derived | Shear modulus for the torsional term, equation (19). |
| `Material.torsion_model` | $J$ | `"rectangle"` | Saint-Venant torsion constant of a solid rectangle, equation (20). |
| `Material.torsion_beta` | $\beta$ | 0.63 | Short-beam torsional stiffening coefficient, equation (21). Zero disables the correction. |
| `Material.resistivity` | $\rho_e$ | none | Electrical resistivity. Sets resistance of beam links, equation (41); without it that resistance is unknown. Only the device layer's value is used. |
| `Material.relative_permittivity` | $\varepsilon_\text{ox}$ | oxide 3.9, silicon none | Static dielectric constant. Only the box layer's value is used, for capacitance through buried oxide, equation (43): silicon is a conductor in the electrical model. |

`Material.silicon()` and `Material.oxide()` return these defaults and accept
keyword overrides.

## `SolveOptions`

| Field | Default | Effect |
| --- | --- | --- |
| `n_modes` | all | Number of lowest modes to return. |
| `method` | `"dense"` | `"dense"` solves the full generalized eigenproblem; `"sparse"` uses shift-invert iteration. |
| `sigma` | 0 | Shift for the sparse solver, in squared angular frequency. |

## `ElectrostaticParams`

| Field | Symbol | Default | Effect |
| --- | --- | --- | --- |
| `ambient_permittivity` | $\varepsilon_a$ | 1 | Relative permittivity of the medium that fills released gaps and trenches, equations (43) and (61). |
| `element_size` | $h$ | none | Longest boundary element on conductor outlines, equation (73). Required when a pair of nets is solved by the field method, in one dense system or in windows; a structure whose pairs are all taken as plates needs none. Half of the narrowest trench is a starting value: computed capacitance rises toward the exact one as it shrinks. |
| `corner_grading` | $g$ | 16 | Ratio of `element_size` to the element length asked for at a convex corner of a right angle or sharper, equation (72). 1 disables refinement. |
| `max_elements` | | 8000 | Most boundary elements solved as one dense system. Memory grows with its square and time with its cube. Pairs asked to be dense whose problem needs more are solved in windows, equations (91) to (96), and reported as windowed; no window grows past this number. |
| `window_elements` | $n_w$ | 500 | Most boundary elements in the core of a window. Each window is one dense solution of its core and the elements around it, and supplies the charge on its core, equation (92). Smaller cores mean more and smaller solutions. |
| `window_tolerance` | $\tau$ | 0.001 | The overlap of a window around its core doubles until the charge on the core changes by less than this fraction, equation (96), or the overlap reaches `interaction_distance`. |
| `dense_cutoff` | $\lambda_d$ | 10 | Pairs of nets whose feature size, the length of their facing walls in units of the gap between them, equation (105), is below this are solved by the field method in one dense system, equation (106). Below ten gaps fringe field carries from an eighth to most of the capacitance. |
| `plates_cutoff` | $\lambda_p$ | 100 | Pairs whose feature size is at least this take their facing walls as parallel plates, equations (97) to (103); pairs between the two cutoffs are solved by the field method in windows. At or above a hundred gaps plates are within about 3.5 percent in capacitance and 1 percent in pull and side spring; see [accuracy and defaults](theory.md#accuracy-and-defaults). Lower it to trade accuracy for time: 15 takes as plates combs whose fingers overlap more than about fifteen gaps. 0 for both cutoffs puts every pair on plates, infinity for both every pair in one dense system. |
| `interaction_distance` | $d_i$ | none | Device electrodes whose silicon lies within this distance, directly or through a chain of neighbors, are taken together for sidewall capacitance, equation (110), and every pair within it is an interacting pair, equation (106). It is also the widest gap across which two walls count as a plate, equation (97), the widest overlap of a window, equation (96), and the reach of the cores of pairs solved by the field. Sidewall capacitance requires it. A few times the widest trench that should count is a starting value. |
| `include_plates` | | true | Whether the capacitance model of a structure holds plates between device and handle electrodes. |
| `include_sidewalls` | | true | Whether it holds sidewall capacitance between device electrodes. Without it neither `element_size` nor `interaction_distance` is needed. |

## `StaticOptions`

| Field | Symbol | Default | Effect |
| --- | --- | --- | --- |
| `tolerance` | $\tau_s$ | 1e-9 | Newton iteration for static equilibrium stops when the strain energy of its last step is below this fraction, squared, of the strain energy of the displacement, equation (122), or one step after it is below this fraction itself. |
| `max_iterations` | | 30 | Most Newton iterations of one solve. Reaching it raises `NoEquilibriumError`. |

## `SmallSignalOptions`

| Field | Symbol | Default | Effect |
| --- | --- | --- | --- |
| `quality_factor` | $Q$ | **required** | Quality factor of every mode of the small-signal model about an equilibrium, equation (132). Damping is an input: the library computes no loss mechanism. It sets the width of each resonance, $f_n / Q$, and the motional resistance of equation (136). |
| `linear_tolerance` | $\tau_l$ | 0.1 | `SmallSignal.linearity` warns with `LinearizationWarning` when the swing across a pair of electrodes against its bias, or the closing of the narrowest gap at the peak of the motion, exceeds this fraction, equation (137). The force then errs by about half the first and the coefficients by about twice the second. |

## `TransientOptions`

| Field | Symbol | Default | Effect |
| --- | --- | --- | --- |
| `quality_factor` | $Q$ | **required** | Quality factor of every mode of the structure without potentials, which sets the damping matrix of a transient, equation (140). Under bias a mode keeps that damping and so has a slightly lower quality factor. |
| `tolerance` | $\tau_t$ | 1e-6 | Error the integrator allows per step, relative to each state, and in absolute terms this times the scales of equation (146): the buried oxide gap for translations, the lowest frequency for velocities, and node capacitance times the largest source value for charges. |
| `max_step` | | none | Longest time step in seconds. Needed only for drive features that error control could step over, such as a short pulse of a callable waveform without breakpoints. |
| `contact_clearance` | $\chi_c$ | 0.1 | The run ends with a contact event when the narrowest gap between two bodies falls to this fraction of its width at rest, equation (126). Nothing in the model describes bodies that touch. |

## `Source`

A voltage source that drives one electrode in a transient.

| Field | Symbol | Default | Effect |
| --- | --- | --- | --- |
| `electrode` | | **required** | Electrode it connects to, as `Structure.electrode_at` returns it for a point on the pad. An integer is refused. |
| `waveform` | $e(t)$ | **required** | Volts against seconds: a number for a constant, a `PiecewiseLinear`, or any continuous callable of time. |
| `resistance` | $R_k$ | 0 | Series resistance in ohms, equation (141). Zero holds the electrode at the waveform. |

## `PiecewiseLinear`

A waveform through points, constant before the first and after the last. Its
points are breakpoints of the integration.

| Field | Symbol | Default | Effect |
| --- | --- | --- | --- |
| `times` | | **required** | Times of the points in seconds, strictly ascending. |
| `values` | | **required** | Value at each time, in volts for a source. |

## `Wire`

A connection made off the die between two pads, passed to `extract` in
`wires`.

| Field | Symbol | Default | Effect |
| --- | --- | --- | --- |
| `start` | | **required** | One end: a point `(x, y)` on a pad of the device layer, or `(x, y, layer)` with layer `"device"` or `"handle"`. The conductor and electrode covering it are those the wire joins. |
| `end` | | **required** | The other end, in the same form. |
| `resistance` | | 0 | Resistance of the wire between its ends, its contacts included. Zero is ideal. It is the resistance of the wire's link. |

## Validation

Each object validates itself on construction and raises `ValueError` naming
the field: lengths, areas, thicknesses, densities, moduli, and resistivity are positive;
`beam_width_min < beam_width_max`; `beam_overlap_margin < body_min_width / 2`;
`torsion_beta >= 0`; `via_resistance >= 0`; `relative_permittivity >= 1`;
`ambient_permittivity >= 1`; `element_size > 0` and `interaction_distance > 0`
when given; `corner_grading >= 1`; `max_elements >= 1`; `window_elements >= 1`;
`window_tolerance > 0`; `dense_cutoff >= 0` and `plates_cutoff >= dense_cutoff`;
`tolerance > 0` and `max_iterations >= 1` of `StaticOptions`;
`quality_factor > 0` and `linear_tolerance > 0` of `SmallSignalOptions`;
`quality_factor > 0`, `tolerance > 0`, `max_step > 0` when given, and
`0 < contact_clearance < 1` of `TransientOptions`; `resistance >= 0` of a
`Source`; `PiecewiseLinear` times strictly ascending, as many as values; a `Wire` end is `(x, y)` or `(x, y, layer)` with layer
`"device"` or `"handle"`, and its `resistance >= 0`; every stack layer has a
thickness.

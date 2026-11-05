# Parameters

Every tunable is a named field on one of four parameter objects. There are no
hidden constants. A test checks that this page and the fields agree. Each
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
| `element_size` | $h$ | none | Longest boundary element on conductor outlines, equation (73). Sidewall capacitance requires it; plates do not use it. Half of the narrowest trench is a starting value: computed capacitance rises toward the exact one as it shrinks. |
| `corner_grading` | $g$ | 16 | Ratio of `element_size` to the element length asked for at a convex corner of a right angle or sharper, equation (72). 1 disables refinement. |
| `max_elements` | | 8000 | Most boundary elements solved as one dense system. Memory grows with its square and time with its cube; a larger mesh is an error. |

## Validation

Each object validates itself on construction and raises `ValueError` naming
the field: lengths, areas, thicknesses, densities, moduli, and resistivity are positive;
`beam_width_min < beam_width_max`; `beam_overlap_margin < body_min_width / 2`;
`torsion_beta >= 0`; `via_resistance >= 0`; `relative_permittivity >= 1`;
`ambient_permittivity >= 1`; `element_size > 0` when given; `corner_grading >= 1`;
`max_elements >= 1`; every stack layer has a thickness.

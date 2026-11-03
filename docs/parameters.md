# Parameters

Every tunable is a named field on one of three configuration objects. There
are no hidden constants. A test checks that this page and the fields agree.

Lengths are in meters, areas in square meters, densities in kilograms per
cubic meter, moduli in pascals. Fields marked **required** have no default:
a wrong release distance or body width yields a plausible but wrong
structure, so the library does not guess them.

## `ExtractionParams`

Symbols refer to the [theory of operation](theory.md).

| Field | Symbol | Default | Effect |
| --- | --- | --- | --- |
| `heal_distance` | $r_h$ | 10 nm | Buffer out then in on device and handle layers to close hairline gaps, equation (1). Beam candidates narrower than $r_h$ are dropped as slivers of the body opening. |
| `simplify_tolerance` | $r_h$ | = `heal_distance` | Douglas–Peucker tolerance applied after healing. |
| `release_hole_max_area` | $A_\text{min}$ | **required** | Holes inside device polygons at or below this area are release holes: they enter the release estimate, are removed from device-layer mass, and are filled for connectivity and beam extraction. |
| `release_distance` | $d$ | **required** | Isotropic undercut of the surviving buried oxide, equation (2). |
| `body_min_width` | $w_\text{min}$ | **required** | Features narrower than this become beam candidates, equation (3). Also the upper width bound of a beam. |
| `beam_overlap_margin` | $\epsilon$ | = `heal_distance` | Expansion of beam candidates so they overlap the bodies they touch, equation (4). |
| `beam_width_min` | | 0 (disabled) | Candidates narrower than this are discarded as slivers. |
| `beam_width_max` | | = `body_min_width` | Candidates wider than this are discarded as non-slender. |
| `anchor` | | largest handle polygon | Override by a point inside the intended handle polygon, or by index into the handle polygons ordered by descending area. |
| `buffer_resolution` | | 8 | Segments per quarter circle of the round release undercut, equation (2). Healing and body opening use mitered corners. |
| `precision_grid` | | off | Optional snap grid applied before overlay operations. |

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

`Material.silicon()` and `Material.oxide()` return these defaults and accept
keyword overrides.

## `SolveOptions`

| Field | Default | Effect |
| --- | --- | --- |
| `n_modes` | all | Number of lowest modes to return. |
| `method` | `"dense"` | `"dense"` solves the full generalized eigenproblem; `"sparse"` uses shift-invert iteration. |
| `sigma` | 0 | Shift for the sparse solver, in squared angular frequency. |

## Validation

Each object validates itself on construction and raises `ValueError` naming
the field: lengths, areas, thicknesses, densities, and moduli are positive;
`beam_width_min < beam_width_max`; `beam_overlap_margin < body_min_width / 2`;
`torsion_beta >= 0`; every stack layer has a thickness.

## TOML files

All three objects load from one file with the same field names:

```toml
[stack.device]
thickness = 25e-6
material = "silicon"

[stack.box]
thickness = 2e-6
material = "oxide"

[stack.handle]
thickness = 400e-6
material = "silicon"

[materials.silicon]
youngs_modulus = 130e9

[extraction]
release_hole_max_area = 100e-12
release_distance = 5e-6
body_min_width = 10e-6

[solve]
n_modes = 20
```

```python
stack, params, options = cf.load_config("soi_25um.toml")
```

# Changelog

Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and
versions follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Electrodes and links: `Structure.electrodes` lists silicon of one conductor
  that moves with one body, `Structure.links` lists beams, bridges, and vias
  that conduct between them, and `Structure.electrode_at` finds an electrode
  by point. Extraction reports an `electrodes` stage and the JSON tree carries
  both lists.
- `Material.resistivity`, `ExtractionParams.via_resistance`, and
  `Link.resistance`: a beam link has resistance from device-layer resistivity
  and the beam's length, width, and thickness; each via link between a device
  and a handle electrode has the supplied via resistance, 0 ohms by default;
  bridges are ideal.
- `Material.relative_permittivity`, 3.9 on the oxide preset and unset on
  silicon.
- `Structure.stack`, `Structure.params`, and `Structure.openings`: the layer
  stack and parameters of the extraction and the release openings that were
  filled in device polygons or drawn on the release layer.
- `Structure.plates` and `Plate`: capacitance between each device electrode
  and the handle electrode beneath it, through surviving oxide and across the
  released gap with release openings removed. `Plate.evaluate` returns
  capacitance, gradient, and Hessian in the degrees of freedom of the two
  bodies: exact in vertical travel, of second order in tilt, and of first
  order in in-plane travel where the handle outline passes under the plate.
  `ElectrostaticParams.ambient_permittivity` sets the medium in the gap.
- `sidewall_capacitance` and `SidewallCapacitance`: capacitance matrix per
  unit thickness between conductors drawn in the layout plane, from a
  boundary element solution of the two-dimensional field with closed-form
  element integrals, with its exact gradient and Hessian for travel and
  rotation of rigid groups of conductors. `ElectrostaticParams` gains `element_size`,
  `corner_grading`, and `max_elements`.
- `Structure.capacitance` and `CapacitanceModel`: capacitance matrix between
  all electrodes at rest, summed from plates and from sidewalls of device
  electrodes solved in clusters, with its oxide, released, and sidewall
  parts and its reduction to nets. `ElectrostaticParams` gains
  `interaction_distance`, `include_plates`, `include_sidewalls`, and
  `skip_large_clusters`.
- `CapacitanceModel.gradient`, `force`, and `stiffness`: change of the
  capacitance matrix in the degrees of freedom of every movable body, and
  the generalized forces and electrostatic stiffness for given electrode
  potentials.

### Removed

- `load_config` and TOML process files. Every parameter object is built in
  Python and passed to the function that uses it, so each analysis takes
  only its own parameters; an application that wants settings files parses
  them itself, for example `ExtractionParams(**table)`.
- `Structure.material`: it duplicated `Structure.stack.device.material`.

## [0.1.0] - 2026-11-04

### Added

- Packaging, code standards, CI, release workflow, and documentation.
- `Material`, `Layer`, `Stack`, `ExtractionParams`, and `SolveOptions` with
  validation, and `load_config` for TOML process files.
- `Layout` with shapely and GDS readers and an adapter registry; polygons
  within a role are merged only where they overlap or touch and stored in
  canonical vertex order. The GDS reader
  rebuilds holes from cut lines and verifies every repaired area against the
  file.
- `extract` producing a `Structure` of rigid bodies and beams: release holes
  identified before healing and again after it, layer healing, oxide estimate, anchor selection, pruning, body assembly
  through the stack, beam qualification with slivers narrower than
  `heal_distance` dropped, stubs attached to their body, a networkx graph, and
  a JSON tree.
- Mass, center of mass, and inertia tensor of every body and its stubs by
  closed-form edge integrals over the stack, with the 6 by 6 body mass matrix.
  Release openings are removed from device-layer mass by subtracting their own
  moments, evaluated in one vectorized pass.
- Beam stiffness: torsion constant with short-beam correction, local 12 by 12
  Euler-Bernoulli matrix, and transformation onto body degrees of freedom.
- Assembly of 6N by 6N mass and stiffness matrices, modal analysis with
  dense and sparse solvers, `ModalResult` with `.npz` export, and a static solve.
- Electrical nets: connected device and handle polygons as conductors, vias
  joining the conductors they overlap, net tags on beams and body polygons,
  `Structure.nets`, and `find_nets`.
- KLayout and gdsfactory readers behind `klayout` and `gdsfactory` extras,
  with an `all` extra for every reader.
- Progress reporting: every stage of `extract`, `modal`, and `static` logs
  start and finish through the `connflex_core` logger, and an optional
  `progress` callback with `start`, `advance`, and `finish` receives the same
  stages, advancing per item where a loop runs.

[Unreleased]: https://github.com/Cao-RQM-Lab/connflex-core/compare/v0.1.0...main
[0.1.0]: https://github.com/Cao-RQM-Lab/connflex-core/releases/tag/v0.1.0

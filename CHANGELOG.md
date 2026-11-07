# Changelog

Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and
versions follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `workers=` on `Structure.capacitance`, `CapacitanceModel.at`,
  `Structure.equilibrium`, `Structure.transient`, and `sidewall_capacitance`:
  threads that solve windows, cut and solve fringe cross-sections, and
  assemble blocks of rows of one dense system, every processor by default. Sums are taken in a fixed
  order, so the count changes results only within rounding.
- `fast` extra: with numba installed, the pair integrals of the field kernel
  and their derivatives run compiled, one pair at a time and outside the
  interpreter's lock; the numpy code remains the fallback, and
  `connflex_core.electrostatics.pairs.COMPILED` tells which is in use.

### Fixed

- A beam candidate that keeps surviving buried oxide underneath, as one wider
  than twice `release_distance` over the handle does, is rejected with the
  oxide area as its reason instead of qualifying as a free beam of full
  length while its oxide was credited to the anchor.
- Two consecutive outline vertices within `1e-6` of `element_size` of each
  other, as healing can leave them, no longer make a boundary element of that
  length, whose row of the field matrix was noise.

### Changed

- Window planning counts the elements of a trial window through the index of
  element centers rather than by testing every element of the mesh.
- Hessian blocks are mapped onto body degrees of freedom by matrix products
  rather than by an `einsum` over six indices.
- Fringe cuts are keyed on a grid of a quarter of the smaller of device
  thickness and released gap, from a sixteenth; their expansion carries the
  rounding to first order, so fewer distinct cross-sections are solved at a
  second-order cost.

## [0.2.0] - 2026-11-06

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
- `Wire` and `extract(..., wires=)`: connections made off the die, such as
  wire bonds and traces, between pads named by a point on them and,
  optionally, their layer. Nets include them, so bonded pads are one net
  in every analysis, and each wire is a link of kind `"wire"` with its
  resistance. `Nets.wires` records the conductors under each end as
  `Connection`.
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
  `interaction_distance`, `include_plates`, and `include_sidewalls`.
- `CapacitanceModel.gradient`, `force`, and `stiffness`: change of the
  capacitance matrix in the degrees of freedom of every movable body, and
  the generalized forces and electrostatic stiffness for given electrode
  potentials.
- `CapacitanceModel.to_dict`: review report of capacitance by pair of
  electrodes and part, between nets, by cluster, and of silicon left out.
- Facing walls as parallel plates: every wall that faces another conductor
  within `interaction_distance` is a plate with it, with exact gradient and
  Hessian for rigid motion. `SidewallCapacitance.walls` and `Cluster.walls`
  hold the pieces of wall, as `FacingWalls`.
- Field solution in windows: cores of at most `window_elements`, each
  solved densely with the elements around it out to an overlap that doubles
  until the charge on the core settles within `window_tolerance` or reaches
  `interaction_distance`. Gradient and Hessian are exact for the windowed
  matrix. `windows` and `asymmetry` on `SidewallCapacitance` and `Cluster`
  report the number of solutions and their mismatch, so no cluster is too
  large to solve.
- Method for each pair of nets: every pair whose electrodes lie within
  `interaction_distance` gets a feature size, the length of its facing walls
  in units of the gap between them, and takes one dense system below
  `ElectrostaticParams.dense_cutoff` (10), windows below `plates_cutoff`
  (100), and plates from there. Pairs are assembled into one matrix with
  exact derivatives. `SidewallCapacitance.feature_size` and `method` and
  `Cluster.interactions` report the choice; a dense problem above
  `max_elements` is solved in windows and reported so.
- `CapacitanceModel.at`: the capacitance model with bodies displaced from
  rest, with matrices, gradient, forces, and stiffness at that displacement.
  Plates are evaluated in closed form and sidewalls on moved outlines, with
  the clusters, mesh, and windows found at rest. `CapacitanceModel.displacement`
  records where a model was evaluated.
- `Structure.equilibrium`, `Equilibrium`, `StaticOptions`, and `NoEquilibriumError`:
  displacement at which beam forces balance electrostatic forces at one set
  of potentials on nets, by Newton iteration on the capacitance model at
  displaced bodies, from rest or from `start`, an earlier equilibrium or a
  displacement. The result carries tangent stiffness, whether the
  equilibrium is stable, and `margin`, the share of beam stiffness left
  where the field softens most. `charges` names floating nets, which keep a
  fixed charge: their potentials follow the bodies, and tangent stiffness
  gains the part that adds. Ramps, pull-in searches, and phase sequences
  are left to the caller, as a series of solves each started from the last.
- `CapacitanceModel.clearance` and `closest`: the narrowest gap between two
  bodies at a displacement as a fraction of its width at rest, and the two
  electrodes across it. `Plate.extent` and `Plate.lowest` give the gap of a
  plate at its closest corner.
- `Structure.small_signal`, `SmallSignal`, and `SmallSignalOptions`: linear
  model about a stable equilibrium for small potentials at ports: `Electrode`
  objects, held by sources, and never integers, which nets also are. Resonant frequencies and mode shapes under
  bias, transduction of each mode at each node, node capacitance with
  floating nets eliminated, and conductance of links with resistance, so
  electrodes behind resistive beams are inner nodes. `admittance` between
  ports and `displacement` for a drive at frequencies the caller gives,
  and `motional` series branches per mode and node. One quality factor for
  every mode is required. `linearity` reports how far a drive swings against the bias and
  closes the narrowest gap, as `Linearity`, and warns with
  `LinearizationWarning` beyond `linear_tolerance`.
- `Equilibrium.floating`: nets that held a fixed charge.
- `SmallSignal.netlist`: SPICE subcircuit of the small-signal model with one
  pin per port, capacitors and resistors between nodes, and one series
  branch per mode coupled to the nodes through controlled sources; all
  modes unless `modes` lists some.
- `Structure.transient`, `Transient`, `TransientOptions`, `Source`, and
  `PiecewiseLinear`: motion, potentials, and charges over time under voltage
  sources at electrodes, each ideal or behind a series resistance, with
  floating nets that keep their charge and the resistor network of links.
  Charges of nodes behind resistance are states beside displacement and
  velocity, integrated by an L-stable implicit method with its Jacobian in
  closed form and restarted at waveform breakpoints. One quality factor
  damps every mode. A run ends at its last time or when a gap closes to
  `contact_clearance`, reported in `contact`. `Energy` is the ledger of
  stored, damped, resistive, and supplied energy at each output time.
- `CurrentSource`, `Resistor`, `Capacitor`, `Switch`, and
  `Structure.transient(..., circuit=)`: current sources into electrodes,
  resistors and capacitors between electrodes or to ground, and switches that
  connect and disconnect a source with resistance or a resistor at given
  times. Nets they reach join the network, and a node that no resistance ties
  to a source or ground starts without charge. The energy ledger counts the
  circuit's capacitors, resistors, and current sources.
- `ElectrostaticParams.surrogate_tolerance` and `Cluster.anchors`: clusters
  solved by the field method reuse their nearest earlier solution at
  displaced bodies, expanded to second order in the motion of their bodies,
  while the expansion stays within the tolerance; the radius around each
  kept solution follows from how well it predicted the next.
  `CapacitanceModel.at(..., exact=True)` solves anew, and an equilibrium
  solves its result anew.
- `ElectrostaticParams.periodic_min_cells`: rows of identical cells in a
  field problem solved in windows, found where outlines repeat their edges,
  take their interior from one cell solved with the periodic kernel and the
  derivatives in turn from one window around it; each cell follows its
  bodies' turn through its lever arm. What a row leaves to the windows is
  cut into cores in the frame of the row, so congruent rows leave windows
  that are translates. On by default for rows of 8 cells or more.
- `ElectrostaticParams.include_fringe`, on by default, `CapacitanceModel.fringe`,
  and `Cluster.sections`: field above and below the device layer along every
  wall, from cuts across walls solved in the plane of their normal and the
  vertical with the handle beneath. It adds to capacitance between device
  electrodes and to the handle, depends on rise and tilt as well as travel
  and turn, and so gives vertical comb springs and the lift of fingers over
  the handle. Facing walls weight it by their overlap, rounded where ends pass
  so that teeth sliding past each other pull without jumps; repeated cuts are
  solved once per structure, and the surrogate expands the fringe with its
  solutions.
- Windows that are translates of one another, element for element with
  their cores, conductors, and bodies, are solved once; each copy moves the
  derivatives to the pivots of its own bodies. Window cores are cut on
  coordinates rounded as for those keys, so translated groups are cut alike.
- Windows and rows of identical cells follow travel: at displaced bodies a
  window keeps its elements while they hold everything within three
  quarters of its overlap of its core, and otherwise takes the elements
  within its overlap where they are now; a cell of a row counts while its
  elements stay inside the runs of the row less `interaction_distance`, and
  the others go to windows. Teeth slid past each other by many pitches keep
  the error of a solution planned where they are.

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

[Unreleased]: https://github.com/Cao-RQM-Lab/connflex-core/compare/v0.2.0...main
[0.2.0]: https://github.com/Cao-RQM-Lab/connflex-core/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/Cao-RQM-Lab/connflex-core/releases/tag/v0.1.0

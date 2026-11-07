# Assumptions and limits

What the model takes for granted, where each assumption enters, what it
costs, and how to check it. Equation numbers refer to the
[theory of operation](theory.md); parameters are those of the
[parameters page](parameters.md). Each entry gives the assumption, its
consequence, and a way to check or mitigate it.

## Layout and extraction

- **Polygons are the fabricated silicon.** No etch bias, no sidewall angle,
  no lag, notching, or scalloping: every layer is a vertical extrusion of its
  polygons, (5) and (157). Beam stiffness in the plane goes as $w^3$, so a
  bias $\delta w$ on a beam of width $w$ errs it by about $3 \delta w / w$,
  and every gap errs by $2 \delta w$. Buffer the polygons by the known bias
  before building the `Layout`; a sloped sidewall cannot be represented.
- **Healing and simplification move outlines by up to $r_h$**, (1): closing
  fills gaps narrower than $2 r_h$, and the Douglas–Peucker pass removes
  features up to `simplify_tolerance` deep and pulls drawn arcs inward by up
  to it. Negligible at the 10 nm default unless gaps are a few hundred
  nanometers; lower `simplify_tolerance` then.
- **Release holes are perforations**, at or below `release_hole_max_area`:
  filled for connectivity and beam candidates, subtracted from the oxide
  estimate (2) and from device-layer mass. A hole inside a beam candidate
  widens the beam it becomes, and a release polygon never severs a beam or a
  body.
- **The undercut is isotropic**, `release_distance` from every oxide edge,
  (2), with a sharp boundary. Real HF undercut depends on time and access;
  the undercut sets which features are released, the oxide under pads, and
  the released part of each plate, (42). Vary `release_distance` and watch
  `Body.box` and `Plate.oxide`.
- **A body is rigid wherever it is at least `body_min_width` wide**, (3),
  and bonded wherever any oxide survives under it. Plate and truss
  compliance is left out; a body held by a sliver of oxide is fully
  clamped. Acute corners sharper than about 23 degrees are beveled by the
  mitre limit of the opening. Watch `Body.box` for slivers.
- **Beams are straight prisms**, (4): a candidate that touches two bodies is
  one Euler-Bernoulli beam of length the chord between its attachments and
  width its area over that length, whatever its shape. A bent, forked, or
  stepped candidate becomes a far stiffer straight beam without warning;
  draw each bend as a block wider than `body_min_width`. A candidate that
  keeps surviving oxide underneath is rejected.
- **Beams carry no mass**; stubs carry mass and no stiffness; candidates
  touching no body or three or more carry neither. Resonances of a mass on
  heavy flexures are high by about half the ratio of flexure mass to proof
  mass. Compare `Body.inertial.mass` with the mass of the flexures.
- **Left out:** residual stress and stress gradient, stiction, root
  compliance and fillets at attachments, aspect-ratio-dependent etch depth,
  non-uniform undercut.

## Mechanics

- **Euler-Bernoulli bending, no shear deformation**, (16) to (19): out of
  plane a beam is too stiff by $1 + 3.07\, (t / L)^2$ for silicon, 5 percent
  at eight thicknesses and 19 percent at four; in the plane the same with
  $w$. Keep flexures slender or expect that error.
- **One isotropic $E$ and $G$ per layout.** Silicon's $E$ runs from 130 to
  169 GPa with orientation and its shear moduli from 51 to 80 GPa; in-plane
  bending of deep beams approaches $E / (1 - \nu^2)$. Set
  `Material.youngs_modulus` and `shear_modulus` for the dominant orientation.
- **The short-beam torsion factor $\beta$ of (21) is unverified**: it has no
  retained source, and its form does not follow restrained-warping theory.
  Set `torsion_beta=0` for the Saint-Venant constant alone when torsion
  matters.
- **Linear small-deflection stiffness at every displacement**, (35): no
  stress stiffening, no geometric nonlinearity, no rotation of the beam frame
  with the bodies, no residual, thermal, or package stress. Results drift from
  reality once a deflection is no longer small against the beam width, or a
  turn against a radian; compare the travel of a result with its flexures.
- **Welded ends on rigid bodies, an immovable anchor**, (27) to (37): no
  root or anchor compliance, no substrate or package motion.
- **Linearized rigid-body kinematics**: rotations are vectors and the mass
  matrix (14) is constant; no gyroscopic or centrifugal terms. Error of order
  $\theta^2$ in the kinematics.
- **Damping is one quality factor for every mode**, `quality_factor` of
  (132) and (140), supplied by the caller: no squeeze-film, slide-film,
  thermoelastic, or anchor loss is computed, and out-of-plane modes in air
  cannot share a $Q$ with in-plane ones. Under bias, (140) keeps the unbiased
  mode shapes, so the damping of a transient and of the small-signal model
  agree only without bias.
- **Loads are six generalized forces at a center of mass**, constant in time;
  no distributed loads, no base acceleration, gravity only as a load.
- **Left out:** thermal expansion and the temperature dependence of $E$,
  fatigue, contact mechanics, stoppers, a nonlinear beam solver.

## Layout-plane field

- **Two-dimensional field times device thickness**, (61) and (111): sidewalls
  are vertical and uniform through the thickness, and the field leaving the
  top and bottom faces enters only through the fringe of (157) to (163).
  Error grows with the gap against the thickness; the fringe table on the
  theory page gives its size.
- **Quasi-static electrostatics in a homogeneous medium**,
  `ambient_permittivity`: perfect conductors, no resistive drop inside
  silicon, no dielectric in trenches, no dielectric, substrate, or surface
  charge. Charging drift of pull-in and offset is outside the model.
- **No ground at infinity**, (63): the conductor set is neutral and the far
  potential floats, so no conductor has capacitance to the package, and the
  coupling is truncated by `interaction_distance` into clusters and windows.
  Walls facing open space are where that truncation costs most.
- **Handle sidewalls, beams, and bridges carry no charge**;
  `CapacitanceModel.unassigned_area` reports the silicon left out.
- **Electrodes move rigidly**, (83) and (118), acting on the field through the
  lever of (114) at the device mid-plane, which is first order in tilt; a
  finger does not bend under its side pull.
- **Piecewise-constant Galerkin elements**, (64) to (74): the energy of one
  dense system is a lower bound and its capacitance rises toward the exact
  one as `element_size` shrinks, at about order 1.2 because of corners. There
  is no error estimate; halve `element_size` and compare. A sum of windows is
  not bounded so.
- **The pair integral is exact to about $10^{-12} L_i L_j$**, (67) to (69),
  with the series beyond eight mean lengths.
- **Windows truncate the far field**, (91) to (96): each core takes the charge
  from a window whose overlap doubles from $\min(4 h, d_i)$ until
  `window_tolerance` is met or `interaction_distance` is reached, so a
  windowed pair whose facing walls are long against $d_i$ keeps an error that
  refinement does not remove, a few percent of its force; `asymmetry` is the
  mismatch of overlapping windows, not that error. Raise
  `interaction_distance` against the facing length, or lower `dense_cutoff`.
- **Window sets and row cells are held while travel allows**, (164) and (165):
  $\mathbf{C}(\mathbf{q})$ steps by about the asymmetry where a set changes, so
  two histories reaching one displacement can differ by one such step.
- **The plates method has no fringe and a cusp where teeth stand tip to tip**,
  (97) to (103): capacitance is low by the fringe share, there is no detent
  stiffness, and (97) jumps where a wall's slope seen from another crosses
  $1/10$. Pairs below `plates_cutoff` are solved by the field.
- **The method of a pair follows its feature size**, (105) and (106), against
  cutoffs calibrated on combs, walls, and toothed racks; a geometry unlike
  those may sit elsewhere in accuracy.
- **Rows of identical cells**, (151) to (156): the periodic kernel takes
  interior cells of combs whose cells are exact translates on outlines whose
  repeating unit has at most 16 edges; the derivatives in turn come from a
  finite window and hold while the turn times half the row length is small
  against the gap.
- **Threads and compiled kernels** change results only by rounding; the
  compiled kernels run with fast-math flags, so a degenerate input that
  validation does not catch returns garbage rather than NaN.
- **Gradient and Hessian are exact derivatives of the discretized, windowed
  matrix**, (76) to (96) and (119): forces are conservative and balance over
  bodies to rounding; their accuracy is that of the mesh.

## Plates and fringe

- **Silicon is a perfect conductor and the box a fixed dielectric**, (42) and
  (43): no depletion on the oxide capacitance of a pad, no dielectric
  charging. Pad capacitance on a lightly doped handle falls with bias in
  reality.
- **The field in the released gap is uniform and vertical**, (43) to (47); the
  plate itself holds no fringing field and release openings carry none, so
  the released part is a lower bound; with `include_fringe` the field at
  device walls over the handle is added by the open cuts of (162).
- **The released gap is the box thickness everywhere**, a flat handle: no
  cavity SOI, no recess, no curl. Plate capacitance goes as $1 / g_0$.
- **Plate kinematics**, (44) to (58): exact in rise, second order in tilt with
  a fourth-order error for symmetric plates, first order in in-plane travel;
  the contact check closes a plate at the corner of its overlap, (126).
- **Flexures have no plate**: beams and bridges belong to no electrode.
- **The fringe is two two-dimensional fields added**, the layout plane and a
  cut across each wall, (157): corners of outlines add no field, and the
  accuracy is tabulated only for $t / g$ of 2.5 and above. The cut holds the
  handle under its whole length or none, decided at one point, and truncates
  the wall's own block at $d_i$, so for blocks wider than about $t$ the fringe
  to a faced conductor grows a little with `interaction_distance`.
- **Cuts are rounded to a grid**, (160): a quarter of the smallest of
  thickness, box, and narrowest gap, and expanded from the rounded cut. The
  grid costs the capacitance second order in the grid, the force about first
  order, and the stiffness first order; the theory page gives measured
  values.
- **A block in a cut rises as a whole**: the slope of a tilted body within the
  plane of the cut is left out, a first-order term of relative size
  $h w / (2 g_0)$ of the fringe to the handle.
- **The cross-section holds air under every block**, even under pads where
  oxide survives beyond the undercut; valid while `release_distance` exceeds
  a few box gaps.
- **Cross-section mesh and slab**: elements at most $t / 8$ on walls, the box
  height near features; the handle slab reaches $d_i + 4 t$ past the cut and
  the excess still moves by about a percent as it grows further.
- **With `surrogate_tolerance`** a cluster is the second-order expansion of
  the nearest kept solution, (148) to (150), with its Hessian frozen: the
  Hessian errs by order $\sqrt{\tau_c}$, which slows Newton and Radau but
  does not move results; the radius is sized by the error of the sidewall
  and fringe expansions.
- **Left out:** tapered walls, handle trenches, a lid electrode, squeeze-film
  damping in the vertical gap, Casimir and van der Waals forces, charging.

## Equilibrium, small signal, transient

- **Quasi-static electrostatics in the coupled problem**: the capacitance at
  each instant is the static field at the displaced bodies; electrical
  dynamics enter only through lumped $R$ and $C$, (139) to (144).
- **Electrodes are equipotential**; resistance exists only in links, (41),
  vias, and wires; a floating net is one node that keeps its charge exactly,
  (127); no leakage through the box.
- **Newton finds the equilibrium near its start**, (121): which of several it
  finds depends on `start`; sequences are the caller's, and a solve from rest
  far beyond pull-in ends in `NoEquilibriumError` without saying whether a
  fold, a far start, or the plates cusp stopped it. With
  `surrogate_tolerance` the displacement converges on the expanded forces and
  carries an error of order $\tau_c$.
- **Stability and margin use the small-deflection tangent**, (123) and (124);
  loose bodies are kept at rest.
- **Nothing describes bodies that touch**: clearance (126) is the gap at the
  closest corner of a plate's overlap or the exact distance between carried
  outlines of electrodes of one cluster; electrodes farther apart than $d_i$
  at rest are never checked. A transient ends at `contact_clearance` and
  misses a bounce through it within one step.
- **The capacitance model at displaced bodies** keeps the clusters, meshes,
  cores, and methods of rest: sidewalls depend on in-plane motion alone,
  rise and tilt enter through the fringe, and with `include_fringe`
  $\mathbf{C}(\mathbf{q})$ steps where a gap crosses a grid line.
- **Small signal is first order about a stable equilibrium**, (129) to (138);
  `linearity` reports the voltage and gap ratios of (137) and only warns. A
  netlist that drops a mode drops its compliance too.
- **The transient integrator** is Radau IIA of order 5 with the tolerance of
  (146); waveforms must be continuous, a callable one is differentiated
  numerically, and a run starts at zero velocity with current-free potentials
  and cannot be resumed from another run's end. The energy ledger of (145)
  measures integrator error plus the steps of the surrogate and the fringe
  grid.
- **Left out:** time-varying mechanical loads and base acceleration, gas
  damping, thermoelastic and anchor loss, adhesion after contact, charging
  drift.

## Numerics and tolerances

Every fixed constant of the methods is listed with its value on the
[parameters page](parameters.md#constants). Tolerances on coordinates are
$10^{-9}$ of a length or reach, rounding grids are $10^{-10}$ and
$10^{-11}$ m, and the drop of residue vertices is $10^{-12}$ m; none ties the
library to a process.

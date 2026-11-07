# Symbols and parameters

Every symbol below is a documented input. Full definitions, units, defaults,
and validation rules are on the [parameters page](../parameters.md).

| Symbol | Meaning | Parameter |
| --- | --- | --- |
| $r_h$ | healing tolerance and simplification tolerance | `ExtractionParams.heal_distance`, `ExtractionParams.simplify_tolerance` |
| $A_\text{min}$ | area threshold for release holes found inside the device layer | `ExtractionParams.release_hole_max_area` |
| $d$ | release undercut distance | `ExtractionParams.release_distance` |
| $w_\text{min}$ | minimum feature width of a rigid body; also the upper width bound of a beam | `ExtractionParams.body_min_width`, `ExtractionParams.beam_width_max` |
| $\epsilon$ | beam candidate expansion | `ExtractionParams.beam_overlap_margin` |
| $t_l$ | layer thickness | `Layer.thickness` on each `Stack` layer |
| $\rho_l$ | layer density | `Material.density` |
| $\rho_e$ | electrical resistivity of the device layer | `Material.resistivity` |
| $R_v$ | resistance of one via link | `ExtractionParams.via_resistance` |
| $g_0$ | gap between device and handle silicon | `Layer.thickness` of the box layer |
| $\varepsilon_\text{ox}$ | relative permittivity of buried oxide | `Material.relative_permittivity` of the box layer |
| $\varepsilon_a$ | relative permittivity of the medium in released gaps and trenches | `ElectrostaticParams.ambient_permittivity` |
| $h$ | longest boundary element | `ElectrostaticParams.element_size` |
| $g$ | refinement of elements at convex corners | `ElectrostaticParams.corner_grading` |
| $d_i$ | distance within which device electrodes are taken together for sidewalls, the widest gap between walls that count as plates, and the widest overlap of a window | `ElectrostaticParams.interaction_distance` |
| $\lambda_d$ | feature size below which a pair of nets is solved in one dense system | `ElectrostaticParams.dense_cutoff` |
| $\lambda_p$ | feature size from which a pair of nets is taken as plates | `ElectrostaticParams.plates_cutoff` |
| $n_w$ | most elements in the core of a window | `ElectrostaticParams.window_elements` |
| $\tau$ | change of the charge on a core below which the overlap of its window stops growing | `ElectrostaticParams.window_tolerance` |
| $\tau_c$ | error of the second-order expansion of a field solution within which displaced clusters reuse it | `ElectrostaticParams.surrogate_tolerance` |
| $n_p$ | fewest identical cells of a row solved by the periodic kernel | `ElectrostaticParams.periodic_min_cells` |
| $\tau_s$ | relative size of the last Newton step at which static equilibrium is accepted | `StaticOptions.tolerance` |
| $Q$ | quality factor of every mode about an equilibrium | `SmallSignalOptions.quality_factor` |
| $\tau_l$ | voltage or gap ratio of a drive above which the small-signal model warns | `SmallSignalOptions.linear_tolerance` |
| $Q$ | quality factor of every mode without potentials in a transient | `TransientOptions.quality_factor` |
| $\tau_t$ | error per step of the transient integrator, relative to each state and to the scales of (146) | `TransientOptions.tolerance` |
| $\chi_c$ | clearance at which a transient ends with a contact event | `TransientOptions.contact_clearance` |
| $R_k$ | series resistance of a source | `Source.resistance` |
| $E$, $G$ | Young's and shear modulus of the device layer | `Material.youngs_modulus`, `Material.shear_modulus` or `Material.poisson_ratio` |
| $\beta$ | short-beam torsional stiffening coefficient | `Material.torsion_beta` |

## Quantities the code reports

| Symbol | Meaning | Equation | In the code |
| --- | --- | --- | --- |
| $m$, $\mathbf{c}$, $\mathbf{I}_b$, $\mathbf{M}_b$ | mass, center of mass, inertia, and mass matrix of a body | (6), (7), (15), (14) | `Body.inertial` |
| $L$, $w$, $t$ | length, width, thickness of a beam | after (4) | `Beam` |
| $\mathbf{M}^\text{full}$, $\mathbf{K}^\text{full}$, $\omega_n$, $\boldsymbol{\phi}_n$ | assembled matrices and modes | (34), (35), (38) | `ModalResult` |
| $C_\text{ox}$, $\varepsilon A / g_0$ | plate capacitance through oxide and across the gap | (43) | `Plate.oxide`, `Plate.released` |
| $\mathbf{c}$, $\partial \mathbf{c} / \partial p$, $\partial^2 \mathbf{c} / \partial p\, \partial q$ | sidewall capacitance per unit thickness and its derivatives | (70), (81), (84) | `SidewallCapacitance.capacitance`, `gradient`, `hessian` |
| asymmetry | mismatch of overlapping windows | after (94) | `SidewallCapacitance.asymmetry`, `Cluster.asymmetry` |
| $\lambda_{AB}$, $m_{AB}$ | feature size and method of a pair | (105), (106) | `SidewallCapacitance.feature_size`, `method`; `Interaction` |
| $\mathbf{C}$ and its parts, $\mathbf{C}^\text{net}$ | capacitance between electrodes and between nets | (112), (113) | `CapacitanceModel.matrix`, `oxide`, `released`, `sidewall`, `fringe`, `net_matrix` |
| $f_k$, $\mathbf{K}_e$, $\boldsymbol{\Gamma}$ | force, electrostatic stiffness, transduction | (116), (129) | `CapacitanceModel.force`, `stiffness`, `charge_gradient` |
| $\chi$ | clearance | (126) | `CapacitanceModel.clearance`, `closest`; `Transient.clearance` |
| $\mathbf{K}_t$, $\lambda_\text{min}$, $\mu_\text{min}$ | tangent stiffness, softest motion, margin | (128), (123), (124) | `Equilibrium.stiffness`, `softest`, `margin`, `stable` |
| $\boldsymbol{\eta}_n$, $\mathbf{C}'$, $\mathbf{G}$, $\mathbf{Y}$ | transduction per mode, node capacitance and conductance, admittance | (132), (131), (130), (135) | `SmallSignal.transduction`, `capacitance`, `conductance`, `admittance` |
| $L_{ns}$, $C_{ns}$, $R_{ns}$ | motional branch | (136) | `SmallSignal.motional` |
| $r_v$, $r_g$ | range of the linearization | (137) | `SmallSignal.linearity` |
| $W$, $W_D$, $W_R$, $W_S$ | energy ledger | (145) | `Transient.energy` |
| $i_h$, $i_k$ | source currents | (144), (141) | `Transient.current` |

# connflex-core

Headless Python library that turns a GDS layout, or any set of polygons, of a
silicon-on-insulator (SOI) MEMS process into a reduced-order electromechanical
model: rigid bodies joined by Euler-Bernoulli beams, their mass and stiffness
matrices, and the resonant modes of the assembled system; the layout's
electrical nets and electrodes; and the capacitance between electrodes with
its derivatives in the motion of every body. From these it solves static
equilibrium under voltage, the small-signal model about it with a SPICE
netlist, and transients under voltage and current sources with a circuit
around the structure.

| Step | Input | Output |
| --- | --- | --- |
| Heal and simplify | device, handle polygons | clean layers |
| Release estimate | device, handle, release holes, undercut distance | surviving oxide regions |
| Body extraction | opened device layer, oxide, handle | rigid bodies anchored at the largest handle polygon |
| Beam extraction | thin device regions bridging two bodies | beams with length, width, thickness |
| Mass properties | body and stub polygons and layer stack | mass, center of mass, inertia tensor |
| Beam stiffness | beam geometry and material | 12 by 12 stiffness on body degrees of freedom |
| Assembly and solve | all bodies and beams | frequencies and mode shapes |
| Electrical nets | device, handle, vias, wires off the die | net graph |
| Electrodes and links | bodies, stubs, beams, vias | nodes and branches of the electrical model |
| Plates to the handle | electrodes, oxide, release openings, stack | capacitance of each device electrode to the handle beneath it, with gradient and Hessian in body degrees of freedom |
| Sidewall capacitance | conductor outlines in the layout plane | capacitance matrix between conductors per unit thickness |
| Capacitance model | electrodes, plates, sidewall clusters | capacitance matrix between electrodes and between nets |
| Fringe | walls of the device layer, its thickness, handle beneath | field above and below the device layer, added between device electrodes and to the handle |
| Static equilibrium | beam stiffness, capacitance model, potentials on nets | displacement under voltage and its stability |
| Small signal | equilibrium, ports, quality factor | modes under bias, admittance, motional branches, SPICE netlist |
| Transient | capacitance model, voltage and current sources, circuit, quality factor | motion, potentials, charges, currents, and energy over time |

Mechanical and electrostatic analyses each run on their own from one extracted
structure. Static equilibrium, small signal, and transients take both.

- [Installation](installation.md)
- [Quick start](quickstart.md)
- [Usage](usage.md)
- [Input conventions](conventions.md)
- [Parameters](parameters.md)
- [Theory of operation](theory.md)
- [Assumptions and limits](assumptions.md)
- [API reference](api/index.md)

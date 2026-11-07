# Theory of operation

These pages are the canonical description of what connflex-core computes: the
model behind every number the library returns, in the order the library builds
it. Each page states what a stage takes and what it produces, gives its
equations, and ends with a table of the symbols it uses against the names they
carry in the code, so that a parameter, a result attribute, or a module can be
found from the mathematics and the other way round.

## How to read these pages

A silicon-on-insulator layout enters as polygons and leaves as a reduced-order
model: rigid bodies joined by straight beams, their mass and stiffness, the
electrodes they carry, and the capacitance between those electrodes as a
function of where the bodies are. The coupled analyses then work on that model
alone. The pages follow that order.

| Page | What it computes | Equations | Entry point in the code |
| --- | --- | --- | --- |
| [Layout, release, and beams](extraction.md) | healed layers, surviving oxide, rigid bodies, beams and rejected candidates | (1) to (4) | `extract` |
| [Mass, stiffness, and modes](mechanics.md) | mass matrix of each body, effective stiffness of each beam, assembled system, resonant modes | (5) to (38) | `extract`, `Structure.modal`, `Structure.static` |
| [Electrodes and links](electrical.md) | nets, electrodes, and the links that conduct between them | (39) to (41) | `extract`, `Structure.nets`, `Structure.electrodes` |
| [Plates between device and handle](plates.md) | capacitance of device silicon to the handle beneath, with its derivatives in the motion of both bodies | (42) to (60) | `Structure.plates`, `Plate.evaluate` |
| [Sidewall capacitance](sidewalls.md) | field between sidewalls in the layout plane, in one dense system or in windows, with derivatives | (61) to (96) | `sidewall_capacitance` |
| [Facing plates and the method for each pair](methods.md) | sidewalls as parallel plates, the feature size of each pair of nets, and which method it takes | (97) to (108) | `sidewall_capacitance`, `Structure.capacitance` |
| [Capacitance of a structure](capacitance.md) | the capacitance matrix between all electrodes, its forces and stiffness, and its evaluation at displaced bodies | (109) to (119) | `Structure.capacitance`, `CapacitanceModel.at` |
| [Static equilibrium](equilibrium.md) | where the structure rests under potentials, its stability and margin, contact, floating nets | (120) to (128) | `Structure.equilibrium` |
| [Small signal](small-signal.md) | linear model about an equilibrium: modes under bias, admittance, motional branches, netlist | (129) to (138) | `Structure.small_signal` |
| [Transient](transient.md) | motion, potentials, and charges in time under sources and a circuit | (139) to (147) | `Structure.transient` |
| [Field solutions over travel](motion.md) | reuse of field solutions at nearby displacements, and rows of identical cells solved as one cell | (148) to (156) | `CapacitanceModel.at` |
| [Fringe above and below the device layer](fringe.md) | field that leaves the top and bottom of every sidewall, across gaps and to the handle | (157) to (163) | `Structure.capacitance` |
| [Windows and rows over travel](travel.md) | how windows and rows follow bodies that travel far | (164), (165) | `CapacitanceModel.at` |
| [Symbols and parameters](symbols.md) | every symbol of these pages against its parameter or result in the code | | |

Equations are numbered in one sequence across the pages, (1) to (165), so that
a number names one equation wherever it is cited: on the
[parameters](../parameters.md) and [assumptions](../assumptions.md) pages, in
the API reference, and in the code. The number range of each page is in the
table above, and a citation of an equation on another page names that page.

## Conventions

Lengths are in meters, masses in kilograms, potentials in volts, capacitance
in farads, and angles in radians; frequencies are reported in hertz. Vectors
and matrices are bold, $\mathbf{q}$ and $\mathbf{K}$; a point of the layout
plane is also written as the complex number $z = x + \mathrm{i} y$ where that
shortens a formula. The top of the device layer is $z = 0$ and the stack
descends from it: device layer, buried oxide, handle. A body has six degrees
of freedom about its center of mass, three translations and three rotations,
in the order $(u_x, u_y, u_z, \theta_x, \theta_y, \theta_z)$, which is the
order of every assembled vector and matrix in the code.

Each page sets apart, in notes titled by their subject, the choices made in
the implementation that the equations themselves do not fix: how a quantity
is evaluated without loss of digits, what is cached, where a tolerance enters.
A reader after the model alone can pass over them.

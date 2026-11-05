# connflex-core

Headless Python library that turns a GDS layout, or any set of polygons, of a
silicon-on-insulator (SOI) MEMS process into a reduced-order mechanical model:
rigid bodies joined by Euler-Bernoulli beams, their mass and stiffness
matrices, and the resonant modes of the assembled system. It also derives the
layout's electrical nets.

| Step | Input | Output |
| --- | --- | --- |
| Heal and simplify | device, handle polygons | clean layers |
| Release estimate | device, handle, release holes, undercut distance | surviving oxide regions |
| Body extraction | opened device layer, oxide, handle | rigid bodies anchored at the largest handle polygon |
| Beam extraction | thin device regions bridging two bodies | beams with length, width, thickness |
| Mass properties | body and stub polygons and layer stack | mass, center of mass, inertia tensor |
| Beam stiffness | beam geometry and material | 12 by 12 stiffness on body degrees of freedom |
| Assembly and solve | all bodies and beams | frequencies and mode shapes |
| Electrical nets | device, handle, vias | net graph |
| Electrodes and links | bodies, stubs, beams, vias | nodes and branches of the electrical model |
| Plates to the handle | electrodes, oxide, release openings, stack | capacitance of each device electrode to the handle beneath it, with gradient and Hessian in body degrees of freedom |
| Sidewall capacitance | conductor outlines in the layout plane | capacitance matrix between conductors per unit thickness |

Mechanical and electrostatic analyses each run on their own from one extracted
structure.

- [Installation](installation.md)
- [Usage](usage.md)
- [Input conventions](conventions.md)
- [Parameters](parameters.md)
- [Theory of operation](theory.md)
- [API reference](api/index.md)

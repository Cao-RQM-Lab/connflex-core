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
| Mass properties | body polygons and layer stack | mass, center of mass, inertia tensor |
| Beam stiffness | beam geometry and material | 12 by 12 stiffness on body degrees of freedom |
| Assembly and solve | all bodies and beams | frequencies and mode shapes |
| Electrical nets | device, handle, vias | net graph |

- [Installation](installation.md)
- [Usage](usage.md)
- [Input conventions](conventions.md)
- [Parameters](parameters.md)
- [Theory of operation](theory.md)
- [API reference](api/index.md)

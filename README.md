# connflex-core

[![CI](https://github.com/Cao-RQM-Lab/connflex-core/actions/workflows/ci.yml/badge.svg)](https://github.com/Cao-RQM-Lab/connflex-core/actions/workflows/ci.yml)
[![Docs](https://github.com/Cao-RQM-Lab/connflex-core/actions/workflows/docs.yml/badge.svg)](https://cao-rqm-lab.github.io/connflex-core/)
[![PyPI](https://img.shields.io/pypi/v/connflex-core.svg)](https://pypi.org/project/connflex-core/)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/Cao-RQM-Lab/connflex-core/blob/main/LICENSE)

Headless Python library that turns a GDS layout, or any set of polygons, of a
silicon-on-insulator (SOI) MEMS process into a reduced-order mechanical model:
rigid bodies joined by Euler-Bernoulli beams, their mass and stiffness
matrices, and the resonant modes of the assembled system. It also derives the
electrical nets of the layout. There is no GUI.

## Installation

Python 3.11 or newer.

```sh
pip install connflex-core          # core
pip install "connflex-core[gds]"   # with the gdstk reader for GDS files
```

## Usage

```python
import connflex_core as cf

layout = cf.load_gds("resonator.gds", layers=cf.LayerMap(device=(1, 0), handle=(2, 0)))
stack = cf.Stack(
    device=cf.Layer(thickness=25e-6, material=cf.Material.silicon()),
    box=cf.Layer(thickness=2e-6, material=cf.Material.oxide()),
    handle=cf.Layer(thickness=400e-6, material=cf.Material.silicon()),
)
params = cf.ExtractionParams(
    release_hole_max_area=100e-12, release_distance=5e-6, body_min_width=10e-6
)

structure = cf.extract(layout, stack, params)
result = structure.modal()
print(result.frequencies_hz[:5])
```

Inputs are four polygon roles: `device`, `handle`, and optionally `release`
and `vias`. Every polygon marks where the named feature is present. Readers
exist for GDS files and shapely geometries, and other readers plug in through
a small adapter registry. The stack, parameters, and solve options can also be
loaded from one TOML file. All units are SI.

## Documentation

<https://cao-rqm-lab.github.io/connflex-core/>

- [Installation](https://cao-rqm-lab.github.io/connflex-core/installation/)
- [Usage](https://cao-rqm-lab.github.io/connflex-core/usage/)
- [Input conventions](https://cao-rqm-lab.github.io/connflex-core/conventions/)
- [Parameters](https://cao-rqm-lab.github.io/connflex-core/parameters/)
- [Theory of operation](https://cao-rqm-lab.github.io/connflex-core/theory/)
- [API reference](https://cao-rqm-lab.github.io/connflex-core/api/)

## Contributing

See [CONTRIBUTING.md](https://github.com/Cao-RQM-Lab/connflex-core/blob/main/CONTRIBUTING.md). Bugs and feature requests go to
[GitHub Issues](https://github.com/Cao-RQM-Lab/connflex-core/issues),
questions to
[GitHub Discussions](https://github.com/Cao-RQM-Lab/connflex-core/discussions).

## AI usage disclosure

Generative AI tools, including Claude Code (Anthropic) and Codex (openAI), were
used in parts of this software. All final design decisions were made by the
authors, and the authors are responsible for the final code and documentation.

## License

MIT, see [LICENSE](https://github.com/Cao-RQM-Lab/connflex-core/blob/main/LICENSE).

# Installation

Python 3.11 or newer.

```sh
pip install connflex-core
```

Readers for layout files and tools, and compiled field kernels, are optional
extras:

| Extra | Installs | Notes |
| --- | --- | --- |
| `gds` | gdstk | Reads GDSII files. |
| `klayout` | klayout | Reads KLayout layouts and cells. GPL-3 licensed. |
| `gdsfactory` | gdsfactory | Reads gdsfactory components. Includes klayout. |
| `fast` | numba | Compiles the field kernels, a few times faster than the numpy code. |
| `all` | all of the above | |

```sh
pip install "connflex-core[gds]"
pip install "connflex-core[all]"
```

From source:

```sh
git clone https://github.com/Cao-RQM-Lab/connflex-core.git
cd connflex-core
python -m venv venv && source venv/bin/activate
python -m pip install -e ".[all]"
```

Verify:

```sh
python -c "import connflex_core; print(connflex_core.__version__)"
```

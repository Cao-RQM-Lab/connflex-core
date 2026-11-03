# Installation

Python 3.11 or newer.

```sh
pip install connflex-core
```

Readers for layout formats are optional extras:

| Extra | Installs | Notes |
| --- | --- | --- |
| `gds` | gdstk | Reads GDSII files. |
| `gdsfactory` | gdsfactory 9.46 or newer | Accepts gdsfactory components. Large dependency tree. |
| `klayout` | klayout | Accepts klayout cells. The klayout package is GPL-3 licensed and is only imported when installed. |
| `all` | all of the above | |

```sh
pip install "connflex-core[gds]"
```

From source:

```sh
git clone https://github.com/Cao-RQM-Lab/connflex-core.git
cd connflex-core
python -m venv venv && source venv/bin/activate
python -m pip install -e ".[gds]"
```

Verify:

```sh
python -c "import connflex_core; print(connflex_core.__version__)"
```

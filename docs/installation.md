# Installation

Python 3.11 or newer.

```sh
pip install connflex-core
```

The GDS reader is an optional extra:

| Extra | Installs | Notes |
| --- | --- | --- |
| `gds` | gdstk | Reads GDSII files. |

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

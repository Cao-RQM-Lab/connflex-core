"""Layout readers.

An adapter is a callable ``adapter(source, **options) -> Layout``. The built-in adapters are
``shapely``, ``gds``, ``klayout``, and ``gdsfactory``. Other packages add adapters with
``register`` or through the ``connflex_core.adapters`` entry-point group.
"""

from __future__ import annotations

from collections.abc import Callable
from importlib import import_module
from importlib.metadata import entry_points

from connflex_core.layout import Layout

Adapter = Callable[..., Layout]

_BUILTIN = {
    "shapely": ("connflex_core.adapters.shapely_", "from_shapely"),
    "gds": ("connflex_core.adapters.gds", "load_gds"),
    "klayout": ("connflex_core.adapters.klayout_", "from_klayout"),
    "gdsfactory": ("connflex_core.adapters.gdsfactory_", "from_gdsfactory"),
}
_REGISTRY: dict[str, Adapter] = {}
_discovered = False


def register(name: str, adapter: Adapter) -> None:
    """Make ``adapter`` available under ``name``."""
    _REGISTRY[name] = adapter


def available() -> list[str]:
    """Names of all adapters."""
    _discover()
    return sorted(_REGISTRY)


def get(name: str) -> Adapter:
    """Return the adapter called ``name``."""
    _discover()
    if name not in _REGISTRY:
        raise KeyError(f"no adapter named {name!r}; available: {', '.join(sorted(_REGISTRY))}")
    return _REGISTRY[name]


def _discover() -> None:
    global _discovered
    if _discovered:
        return
    _discovered = True
    for name, (module, attribute) in _BUILTIN.items():
        _REGISTRY.setdefault(name, getattr(import_module(module), attribute))
    for entry_point in entry_points(group="connflex_core.adapters"):
        _REGISTRY.setdefault(entry_point.name, entry_point.load())

"""Reader for gdsfactory components."""

from __future__ import annotations

from typing import TYPE_CHECKING

from connflex_core.adapters.klayout_ import from_klayout
from connflex_core.layout import LayerMap, Layout

if TYPE_CHECKING:
    import gdsfactory


def from_gdsfactory(component: gdsfactory.Component, *, layers: LayerMap) -> Layout:
    """Read a gdsfactory component into a layout.

    Database units are applied and instances are flattened.
    """
    return from_klayout(component.kdb_cell, layers=layers)

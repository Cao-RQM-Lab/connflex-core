"""Reduced-order mechanical models and electrical nets from SOI layouts."""

from importlib.metadata import version

from connflex_core import adapters
from connflex_core.adapters.gds import load_gds
from connflex_core.config import load_config
from connflex_core.layout import LayerMap, Layout
from connflex_core.materials import Material
from connflex_core.params import ExtractionParams, SolveOptions
from connflex_core.stack import Layer, Stack

__version__ = version("connflex-core")

__all__ = [
    "ExtractionParams",
    "Layer",
    "LayerMap",
    "Layout",
    "Material",
    "SolveOptions",
    "Stack",
    "__version__",
    "adapters",
    "load_config",
    "load_gds",
]

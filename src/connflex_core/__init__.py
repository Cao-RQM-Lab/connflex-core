"""Reduced-order mechanical models and electrical nets from SOI layouts."""

from importlib.metadata import version

from connflex_core.config import load_config
from connflex_core.materials import Material
from connflex_core.params import ExtractionParams, SolveOptions
from connflex_core.stack import Layer, Stack

__version__ = version("connflex-core")

__all__ = [
    "ExtractionParams",
    "Layer",
    "Material",
    "SolveOptions",
    "Stack",
    "__version__",
    "load_config",
]

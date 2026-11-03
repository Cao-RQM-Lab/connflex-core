"""Reduced-order mechanical models and electrical nets from SOI layouts."""

from importlib.metadata import version

from connflex_core.materials import Material
from connflex_core.stack import Layer, Stack

__version__ = version("connflex-core")

__all__ = ["Layer", "Material", "Stack", "__version__"]

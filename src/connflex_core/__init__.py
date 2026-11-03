"""Reduced-order mechanical models and electrical nets from SOI layouts."""

from importlib.metadata import version

from connflex_core import adapters
from connflex_core.adapters.gds import load_gds
from connflex_core.config import load_config
from connflex_core.extract import extract
from connflex_core.layout import LayerMap, Layout
from connflex_core.materials import Material
from connflex_core.model import Beam, Body, RejectedCandidate, Structure
from connflex_core.params import ExtractionParams, SolveOptions
from connflex_core.stack import Layer, Stack

__version__ = version("connflex-core")

__all__ = [
    "Beam",
    "Body",
    "ExtractionParams",
    "Layer",
    "LayerMap",
    "Layout",
    "Material",
    "RejectedCandidate",
    "SolveOptions",
    "Stack",
    "Structure",
    "__version__",
    "adapters",
    "extract",
    "load_config",
    "load_gds",
]

"""Reduced-order electromechanical models of SOI MEMS from layouts."""

from importlib.metadata import version

from connflex_core import adapters
from connflex_core.adapters.gds import load_gds
from connflex_core.adapters.gdsfactory_ import from_gdsfactory
from connflex_core.adapters.klayout_ import from_klayout
from connflex_core.coupled.small_signal import Linearity, LinearizationWarning, SmallSignal
from connflex_core.coupled.static import Equilibrium, NoEquilibriumError
from connflex_core.coupled.transient import Energy, Transient
from connflex_core.electrical.drivers import (
    Capacitor,
    CurrentSource,
    PiecewiseLinear,
    Resistor,
    Source,
    Switch,
)
from connflex_core.electrical.electrodes import Electrode, Link
from connflex_core.electrical.nets import Conductor, Connection, Nets, Via, Wire, find_nets
from connflex_core.electrostatics.capacitance import (
    CapacitanceModel,
    Cluster,
    HessianBlock,
    Interaction,
)
from connflex_core.electrostatics.facing import FacingWalls
from connflex_core.electrostatics.plates import Plate
from connflex_core.electrostatics.sidewalls import SidewallCapacitance, sidewall_capacitance
from connflex_core.extract import extract
from connflex_core.layout import LayerMap, Layout
from connflex_core.materials import Material
from connflex_core.mechanics.inertia import MassProperties
from connflex_core.mechanics.modal import ModalResult
from connflex_core.model import Beam, Body, RejectedCandidate, Structure
from connflex_core.params import (
    ElectrostaticParams,
    ExtractionParams,
    SmallSignalOptions,
    SolveOptions,
    StaticOptions,
    TransientOptions,
)
from connflex_core.progress import Progress
from connflex_core.stack import Layer, Stack

__version__ = version("connflex-core")

__all__ = [
    "Beam",
    "Body",
    "CapacitanceModel",
    "Capacitor",
    "Cluster",
    "Conductor",
    "Connection",
    "CurrentSource",
    "Electrode",
    "ElectrostaticParams",
    "Energy",
    "Equilibrium",
    "ExtractionParams",
    "FacingWalls",
    "HessianBlock",
    "Interaction",
    "Layer",
    "LayerMap",
    "Layout",
    "Linearity",
    "LinearizationWarning",
    "Link",
    "MassProperties",
    "Material",
    "ModalResult",
    "Nets",
    "NoEquilibriumError",
    "PiecewiseLinear",
    "Plate",
    "Progress",
    "RejectedCandidate",
    "Resistor",
    "SidewallCapacitance",
    "SmallSignal",
    "SmallSignalOptions",
    "SolveOptions",
    "Source",
    "Stack",
    "StaticOptions",
    "Structure",
    "Switch",
    "Transient",
    "TransientOptions",
    "Via",
    "Wire",
    "__version__",
    "adapters",
    "extract",
    "find_nets",
    "from_gdsfactory",
    "from_klayout",
    "load_gds",
    "sidewall_capacitance",
]

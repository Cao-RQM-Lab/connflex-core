"""Process description files."""

from __future__ import annotations

import dataclasses
import tomllib
from collections.abc import Callable, Mapping
from os import PathLike
from typing import TYPE_CHECKING, Any, TypeVar

from connflex_core.materials import Material
from connflex_core.params import ExtractionParams, SolveOptions
from connflex_core.stack import Layer, Stack

if TYPE_CHECKING:
    from _typeshed import DataclassInstance

T = TypeVar("T", bound="DataclassInstance")

PRESETS: dict[str, Callable[..., Material]] = {"silicon": Material.silicon, "oxide": Material.oxide}
ROLES = ("device", "box", "handle")


def load_config(path: str | PathLike[str]) -> tuple[Stack, ExtractionParams, SolveOptions]:
    """Load a stack, extraction parameters, and solve options from one TOML file.

    The file holds the tables ``[stack.device]``, ``[stack.box]``, ``[stack.handle]``, and
    ``[extraction]``, and optionally ``[solve]`` and ``[materials.<name>]``. Each stack table
    gives a ``thickness`` and names its ``material``, which is a preset (``silicon`` or
    ``oxide``) or a ``[materials.<name>]`` table. A materials table overrides fields of the
    preset of the same name or defines a new material. Unknown keys are errors.
    """
    with open(path, "rb") as f:
        data = tomllib.load(f)
    _reject_unknown(data, {"stack", "materials", "extraction", "solve"}, "the top level")
    materials = {
        name: _material(name, table)
        for name, table in _table(data, "materials", required=False).items()
    }
    stack_table = _table(data, "stack")
    _reject_unknown(stack_table, set(ROLES), "[stack]")
    layers = {role: _layer(role, _table(stack_table, role, "[stack]"), materials) for role in ROLES}
    stack = Stack(**layers)
    params = _build(ExtractionParams, _table(data, "extraction"), "[extraction]")
    options = _build(SolveOptions, _table(data, "solve", required=False), "[solve]")
    return stack, params, options


def _table(
    data: Mapping[str, Any], key: str, context: str = "the top level", *, required: bool = True
) -> dict[str, Any]:
    if key not in data:
        if required:
            raise ValueError(f"{context} is missing the [{key}] table")
        return {}
    table = data[key]
    if not isinstance(table, dict):
        raise ValueError(f"{key} in {context} must be a table")
    return table


def _reject_unknown(table: Mapping[str, Any], allowed: set[str], context: str) -> None:
    unknown = sorted(set(table) - allowed)
    if unknown:
        raise ValueError(f"unknown key(s) in {context}: {', '.join(unknown)}")


def _build(cls: type[T], table: Mapping[str, Any], context: str) -> T:
    fields = dataclasses.fields(cls)
    _reject_unknown(table, {f.name for f in fields}, context)
    missing = [
        f.name
        for f in fields
        if f.default is dataclasses.MISSING
        and f.default_factory is dataclasses.MISSING
        and f.name not in table
    ]
    if missing:
        raise ValueError(f"{context} is missing {', '.join(missing)}")
    values = {k: tuple(v) if isinstance(v, list) else v for k, v in table.items()}
    make: Callable[..., T] = cls
    return make(**values)


def _material(name: str, table: Mapping[str, Any]) -> Material:
    context = f"[materials.{name}]"
    if name in PRESETS:
        _reject_unknown(table, {f.name for f in dataclasses.fields(Material)}, context)
        return PRESETS[name](**table)
    return _build(Material, table, context)


def _layer(role: str, table: Mapping[str, Any], materials: Mapping[str, Material]) -> Layer:
    context = f"[stack.{role}]"
    _reject_unknown(table, {"thickness", "material"}, context)
    if "thickness" not in table:
        raise ValueError(f"{context} is missing thickness")
    name = table.get("material")
    if not isinstance(name, str):
        raise ValueError(f"{context} must name its material")
    if name in materials:
        material = materials[name]
    elif name in PRESETS:
        material = PRESETS[name]()
    else:
        raise ValueError(f"unknown material {name!r} in {context}")
    return Layer(thickness=table["thickness"], material=material)

"""Field of a mesh too large for one dense solution, solved in overlapping windows."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import shapely
from numpy.typing import NDArray
from shapely import STRtree

from connflex_core.electrostatics.dense import Charges, charge_derivatives, solve_charges
from connflex_core.electrostatics.outlines import Mesh
from connflex_core.progress import quiet

Array = NDArray[np.float64]
Index = NDArray[np.intp]

FIRST_OVERLAP = 4.0


def cut_cores(mesh: Mesh, most: int) -> list[Index]:
    """Split the elements of a mesh into compact groups of at most ``most``, equation (91).

    A group is halved at the median of its element centers along the longer side of
    their bounding box until it is small enough.
    """
    center = (mesh.start + mesh.end) / 2
    cores, pending = [], [np.arange(len(mesh.conductor))]
    while pending:
        group = pending.pop()
        if len(group) <= most:
            cores.append(group)
            continue
        points = center[group]
        axis = int(np.argmax(np.ptp(points, axis=0)))
        order = group[np.argsort(points[:, axis], kind="stable")]
        pending += [order[len(group) // 2 :], order[: len(group) // 2]]
    return cores


def solve_windows(
    mesh: Mesh,
    body: Index,
    pivots: Array,
    cores: list[Index],
    size: float,
    reach: float | None,
    tolerance: float,
    limit: int,
    advance: Callable[[], None],
) -> tuple[Array, Array, Array, Array]:
    """Sum over windows of the charge that each solution puts on its core.

    A window holds a core and every element within an overlap of its bounding box,
    equations (91) to (96). The overlap starts at ``FIRST_OVERLAP`` times ``size`` and
    doubles until the charge on the core changes by less than ``tolerance``, the overlap
    reaches ``reach``, or the next window would hold more than ``limit`` elements.
    Returns capacitance, gradient, Hessian, and element charges per unit permittivity,
    in the layout of ``SidewallCapacitance`` and before capacitance is made symmetric.
    ``advance`` is called once per core.
    """
    count, conductors, rigid = len(mesh.conductor), len(body), len(pivots)
    center = (mesh.start + mesh.end) / 2
    tree = STRtree(np.asarray(shapely.points(center)))
    capacitance = np.zeros((conductors, conductors))
    gradient = np.zeros((conductors, conductors, rigid, 3))
    hessian = np.zeros((conductors, conductors, rigid, 3, rigid, 3))
    charge = np.zeros((count, conductors))
    three = np.arange(3)
    for core in cores:
        low, high = center[core].min(axis=0), center[core].max(axis=0)
        overlap = FIRST_OVERLAP * size if reach is None else min(FIRST_OVERLAP * size, reach)
        held: tuple[Index, Index, Charges] | None = None
        share = np.zeros((conductors, conductors))
        smallest = True
        while True:
            inside = np.sort(tree.query(shapely.box(*(low - overlap), *(high + overlap))))
            if len(inside) > limit:
                if smallest:
                    raise ValueError(
                        f"a window needs {len(inside)} boundary elements, more than "
                        f"max_elements = {limit}: raise max_elements or element_size, or "
                        f"lower window_elements"
                    )
                break
            smallest = False
            present = np.unique(mesh.conductor[inside])
            before, share = share, np.zeros((conductors, conductors))
            settled = False
            if len(present) > 1:
                part = Mesh(
                    mesh.start[inside],
                    mesh.end[inside],
                    np.searchsorted(present, mesh.conductor[inside]),
                )
                charges = solve_charges(part, len(present), quiet, np.isin(inside, core))
                share[np.ix_(present, present)] = charges.capacitance
                settled = held is not None and bool(
                    np.max(np.abs(share - before)) < tolerance * np.max(np.abs(share))
                )
                held = (inside, present, charges)
            if settled or len(inside) == count or (reach is not None and overlap >= reach):
                break
            overlap = 2 * overlap if reach is None else min(2 * overlap, reach)
        if held is not None:
            inside, present, charges = held
            capacitance[np.ix_(present, present)] += charges.capacitance
            moving = np.unique(body[present])
            slope, bend = charge_derivatives(
                charges, np.searchsorted(moving, body[present]), pivots[moving], quiet
            )
            gradient[np.ix_(present, present, moving)] += slope
            hessian[np.ix_(present, present, moving, three, moving, three)] += bend
            mine = np.isin(inside, core)
            charge[np.ix_(inside[mine], present)] = charges.charge[mine]
        advance()
    return capacitance, gradient, hessian, charge

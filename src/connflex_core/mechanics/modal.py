"""Modal and static solutions."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from os import PathLike

import numpy as np
import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg
from numpy.typing import NDArray

from connflex_core.mechanics.assembly import Dof, assemble
from connflex_core.model import Structure
from connflex_core.params import SolveOptions


@dataclass(frozen=True, eq=False)
class ModalResult:
    """Resonant modes of a structure, equation (38).

    Attributes:
        frequencies_hz: Resonant frequencies in ascending order.
        modes: One column per mode over the rows of ``dofs``, mass-normalized.
        mass: Assembled mass matrix.
        stiffness: Assembled stiffness matrix.
        dofs: (body index, degree of freedom) for each row: ux, uy, uz, tx, ty, tz.
    """

    frequencies_hz: NDArray[np.float64]
    modes: NDArray[np.float64]
    mass: NDArray[np.float64]
    stiffness: NDArray[np.float64]
    dofs: tuple[Dof, ...]

    def to_npz(self, path: str | PathLike[str]) -> None:
        """Save the arrays plus ``dof_body`` and ``dof_name`` to a ``.npz`` file."""
        np.savez(
            path,
            frequencies_hz=self.frequencies_hz,
            modes=self.modes,
            mass=self.mass,
            stiffness=self.stiffness,
            dof_body=np.array([body for body, _ in self.dofs]),
            dof_name=np.array([name for _, name in self.dofs]),
        )


def modal(structure: Structure, options: SolveOptions | None = None) -> ModalResult:
    """Solve the generalized eigenproblem of the assembled matrices for resonant modes.

    The dense solver returns all modes or the lowest ``n_modes``. The sparse solver uses
    shift-invert iteration about ``sigma`` and needs ``n_modes`` smaller than the number of
    degrees of freedom and a nonsingular shifted stiffness.
    """
    options = options or SolveOptions()
    system = assemble(structure)
    size = system.mass.shape[0]
    if size == 0:
        raise ValueError("structure has no movable bodies")
    count = size if options.n_modes is None else min(options.n_modes, size)
    if options.method == "sparse":
        if count >= size:
            raise ValueError("sparse solver needs n_modes below the number of degrees of freedom")
        values, vectors = scipy.sparse.linalg.eigsh(
            scipy.sparse.csr_matrix(system.stiffness),
            k=count,
            M=scipy.sparse.csr_matrix(system.mass),
            sigma=options.sigma,
        )
        order = np.argsort(values)
        values, vectors = values[order], vectors[:, order]
    else:
        values, vectors = scipy.linalg.eigh(
            system.stiffness, system.mass, subset_by_index=[0, count - 1]
        )
    frequencies = np.sqrt(np.clip(values, 0.0, None)) / (2 * np.pi)
    return ModalResult(
        frequencies_hz=frequencies,
        modes=vectors,
        mass=system.mass,
        stiffness=system.stiffness,
        dofs=system.dofs,
    )


def static(structure: Structure, loads: Mapping[int, Sequence[float]]) -> NDArray[np.float64]:
    """Displacements under static loads.

    ``loads`` maps a body index to its six generalized forces (fx, fy, fz, tx, ty, tz) at its
    center of mass. The result has one entry per row of the assembled matrices.
    """
    system = assemble(structure)
    rows = {body: row for row, (body, name) in enumerate(system.dofs) if name == "ux"}
    force = np.zeros(len(system.dofs))
    for body, load in loads.items():
        if body not in rows:
            raise ValueError(f"body {body} is the anchor or does not exist")
        if len(load) != 6:
            raise ValueError(f"load on body {body} must have six components")
        force[rows[body] : rows[body] + 6] = load
    try:
        displacement: NDArray[np.float64] = scipy.linalg.solve(
            system.stiffness, force, assume_a="sym"
        )
    except scipy.linalg.LinAlgError as exc:
        raise ValueError("stiffness matrix is singular: a body has no path to the anchor") from exc
    return displacement

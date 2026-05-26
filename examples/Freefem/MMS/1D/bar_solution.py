"""Discrete 1D FEM solution returned by the solver."""

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class BarSolution1D:
    """Snapshot of a 1D bar FEM solution.

    x0    : rest positions, shape (nx,)
    edges : connectivity,   shape (n_edges, 2)
    u_h   : nodal displacement, shape (nx,)
    """
    x0:    np.ndarray
    edges: np.ndarray
    u_h:   np.ndarray

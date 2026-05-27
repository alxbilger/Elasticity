"""Snapshot of one 2D solve: mesh + connectivity + displacement field."""

from dataclasses import dataclass

import numpy as np


@dataclass
class BeamSolution2D:
    nodes : np.ndarray   # (N, 2) or (N, 3)
    conn  : np.ndarray   # element connectivity (element-type-specific)
    ux    : np.ndarray   # (N,) x-displacement
    uy    : np.ndarray   # (N,) y-displacement

"""Snapshot of one 3D solve: mesh + connectivity + displacement field."""

from dataclasses import dataclass

import numpy as np


@dataclass
class SolidSolution3D:
    nodes : np.ndarray   # (N, 3)
    conn  : np.ndarray   # (n_hex, 8)  hexahedral connectivity from SOFA
    ux    : np.ndarray   # (N,) x-displacement
    uy    : np.ndarray   # (N,) y-displacement
    uz    : np.ndarray   # (N,) z-displacement

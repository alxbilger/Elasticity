"""Abstract base class for 3D MMS cases (Cartesian domain [0,L]^3)."""

import os
import sys
from abc import abstractmethod
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mms_case import MMSCase


def lame(E, nu):
    """Lamé parameters (lambda, mu) for 3D linear elasticity (full Hooke)."""
    lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    mu  = E / (2.0 * (1.0 + nu))
    return lam, mu


class MMSCase3D(MMSCase):
    # Body-force quadrature rule — must be set by each concrete case.
    # No framework fallback; assembly raises if unset.
    source_quadrature_hex = None   # element rule for Q1 hexes (e.g. hex_q1_rule(2))

    @abstractmethod
    def u_ex(self, x, y, z, L):
        """Exact solution: returns (ux, uy, uz)."""

    @abstractmethod
    def grad_u_ex(self, x, y, z, L):
        """Exact 3x3 displacement gradient
        [[dux/dx, dux/dy, dux/dz],
         [duy/dx, duy/dy, duy/dz],
         [duz/dx, duz/dy, duz/dz]]."""

    @abstractmethod
    def source(self, x, y, z, E, nu, L):
        """Body force: returns (fx, fy, fz)."""

    @abstractmethod
    def apply_bcs(self, Solid, nodes_3d, L):
        """Install Dirichlet constraints on the SOFA `Solid` node.

        The MMS author chooses the BC pattern matching the manufactured
        solution. Neumann tractions on the six faces are already assembled
        into the consistent nodal force by the framework (via `self.traction`);
        this method only needs to add SOFA constraint objects.
        """

    def traction(self, x, y, z, nx, ny, nz, E, nu, L):
        """Traction on a face with outward unit normal (nx, ny, nz):
        returns (Tx, Ty, Tz).

        Derived from `grad_u_ex` via sigma = lambda tr(eps) I + 2 mu eps and
        T = sigma . n.
        """
        lam, mu = lame(E, nu)
        G = self.grad_u_ex(x, y, z, L)
        eps = 0.5 * (G + G.T)
        tr  = eps[0, 0] + eps[1, 1] + eps[2, 2]
        S = lam * tr * np.eye(3) + 2.0 * mu * eps
        T = S @ np.array([nx, ny, nz])
        return T[0], T[1], T[2]

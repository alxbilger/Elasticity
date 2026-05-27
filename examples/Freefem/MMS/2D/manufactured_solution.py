"""Abstract base class for 2D MMS cases (Cartesian domain [0,L]^2)."""

from abc import ABC, abstractmethod


def lame(E, nu, dim):
    """Lamé parameters (lambda, mu) for linear elasticity.

    dim="2d" → plane stress:  lambda = E nu / (1 - nu^2)
    dim="3d" → plane strain:  lambda = E nu / ((1 + nu)(1 - 2 nu))
    mu = E / (2 (1 + nu))   (same in both branches)
    """
    if dim == "2d":
        lam = E * nu / (1.0 - nu**2)
    else:
        lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    mu = E / (2.0 * (1.0 + nu))
    return lam, mu


class MMSCase2D(ABC):
    name       = None  # case identifier (must match the params.json key)
    plot_label = None  # LaTeX label for the exact solution

    # Body-force quadrature rules — must be set by each concrete case.
    # No framework fallback; assembly raises if either is unset.
    source_quadrature_quad = None   # element rule for Q1 quads (e.g. quad_q1_rule(2))
    source_quadrature_tri  = None   # element rule for P1 triangles (e.g. tri_p1_rule(3))

    @abstractmethod
    def u_ex(self, x, y, L):
        """Exact solution: returns (ux, uy)."""

    @abstractmethod
    def grad_u_ex(self, x, y, L):
        """Exact displacement gradient: returns 2x2 array
        [[dux/dx, dux/dy], [duy/dx, duy/dy]]."""

    @abstractmethod
    def source(self, x, y, E, nu, L, dim):
        """Body force: returns (fx, fy).

        `dim` selects the constitutive branch:
            "2d" → plane stress, "3d" → plane strain.
        """

    @abstractmethod
    def apply_bcs(self, Beam, nodes_2d, L, dim):
        """Install Dirichlet BCs (and any case-specific constraints) on `Beam`.

        The MMS author chooses the BC pattern matching the manufactured
        solution. Neumann tractions on the four edges are already assembled
        into the consistent nodal force by the framework (via `self.traction`);
        this method only needs to add SOFA constraint objects.
        """

    def traction(self, x, y, nx, ny, E, nu, L, dim):
        """Traction on a face with outward unit normal (nx, ny): returns (Tx, Ty).

        Derived from `grad_u_ex` via sigma = lambda tr(eps) I + 2 mu eps and
        T = sigma . n. `dim` selects plane stress / plane strain (see `source`).
        """
        lam, mu = lame(E, nu, dim)
        G = self.grad_u_ex(x, y, L)
        dux_dx, dux_dy = G[0, 0], G[0, 1]
        duy_dx, duy_dy = G[1, 0], G[1, 1]
        sxx = (lam + 2.0 * mu) * dux_dx + lam * duy_dy
        syy = lam * dux_dx + (lam + 2.0 * mu) * duy_dy
        sxy = mu * (dux_dy + duy_dx)
        return (sxx * nx + sxy * ny,
                sxy * nx + syy * ny)

"""Abstract base class for 1D MMS cases (non-dimensional, x ∈ [0,1])."""

from abc import ABC, abstractmethod


class MMSCase1D(ABC):
    name              = None  # case identifier (must match the params.json key)
    plot_label        = None  # LaTeX label for the exact solution
    source_quadrature = None  # quadrature rule for body-force assembly only

    @abstractmethod
    def u_ex(self, xi):
        """Exact solution."""

    @abstractmethod
    def du_ex(self, xi):
        """Exact first derivative."""

    @abstractmethod
    def source(self, xi, E):
        """Body force f(x)."""

    def traction_bc(self, E):
        """Neumann traction at x=1: F_N = E · u'(1)."""
        return E * self.du_ex(1.0)

    @abstractmethod
    def apply_bcs(self, Bar, E_eff, nx):
        """Install Dirichlet + Neumann BCs on the SOFA `Bar` node.

        The MMS author chooses the BC pattern matching the manufactured
        solution. Typical: Dirichlet at node 0, Neumann tip force at node
        nx-1 using `self.traction_bc(E_eff)`.
        """

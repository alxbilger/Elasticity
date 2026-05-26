"""FEM toolbox: quadrature, assembly, error norms.
"""
import numpy as np

# ---------------------------------------------------------------------------
# Quadrature rules
# Approximation: integral_{x_i}^{x_{i+1}} g(x) dx
#                ~= (h/2) * sum_k( w_k * g(x_k) )
# where x_k = (x_i + x_{i+1})/2 + (h/2)*xi_k
# ---------------------------------------------------------------------------

# Canonical Gauss-Legendre points and weights on the reference segment [-1, 1].
_GAUSS_LEGENDRE_1D = {
    1: (np.array([0.0]),
        np.array([2.0])),
    2: (np.array([-1.0 / np.sqrt(3.0), 1.0 / np.sqrt(3.0)]),
        np.array([1.0, 1.0])),
}


def line_quadrature(n_pts):
    """Return a Gauss-Legendre rule on a physical segment using `n_pts` points.

    The returned callable approximates int_{x1}^{x2} g(x) dx.
    """
    if n_pts not in _GAUSS_LEGENDRE_1D:
        raise ValueError(f"line_quadrature: {n_pts}-point rule not supported")
    xi, w = _GAUSS_LEGENDRE_1D[n_pts]

    def rule(g, x1, x2):
        h   = x2 - x1
        x_k = 0.5 * (x1 + x2) + 0.5 * h * xi
        return 0.5 * h * sum(w_i * g(x_i) for w_i, x_i in zip(w, x_k))
    return rule


# Error-norm rules. H1 requires >=2 pts: the P1 gradient superconverges at
# the element midpoint, so 1-pt fakes O(h^2) instead of O(h).
L2_QUADRATURE = line_quadrature(2)
H1_QUADRATURE = line_quadrature(2)


# ---------------------------------------------------------------------------
# FEM assembly
# ---------------------------------------------------------------------------

def assemble_nodal_forces(f_body, nodes, edges, quadrature):
    """
    Assemble the consistent nodal force vector F_i = integral f_body(x) phi_i(x) dx.

    f_body    : callable x -> float
    nodes     : 1-D array of node coordinates
    edges     : iterable of (a, b) node-index pairs defining each element
    quadrature: rule from line_quadrature(n_pts)
    """
    forces = np.zeros(len(nodes))
    for a, b in edges:
        x1, x2 = nodes[a], nodes[b]
        h = x2 - x1
        forces[a] += quadrature(lambda x, x1=x1, x2=x2, h=h: f_body(x) * (x2 - x) / h, x1, x2)
        forces[b] += quadrature(lambda x, x1=x1, x2=x2, h=h: f_body(x) * (x - x1) / h, x1, x2)
    return forces


# ---------------------------------------------------------------------------
# Error norms
# ---------------------------------------------------------------------------

def l2_error(nodes, edges, u_h, u_ex, quadrature):
    """L2 error norm: sqrt( integral (u_h - u_ex)^2 dx ) over the mesh."""
    total = 0.0
    for a, b in edges:
        x1, x2 = nodes[a], nodes[b]
        h = x2 - x1
        u_a, u_b = u_h[a], u_h[b]
        u_interp = lambda x, x1=x1, h=h, u_a=u_a, u_b=u_b: u_a + (u_b - u_a) * (x - x1) / h
        total += quadrature(lambda x: (u_interp(x) - u_ex(x)) ** 2, x1, x2)
    return np.sqrt(total)


def h1_semi_error(nodes, edges, u_h, du_ex, quadrature):
    """H1 semi-norm error: sqrt( integral (du_h - du_ex)^2 dx ) over the mesh."""
    total = 0.0
    for a, b in edges:
        x1, x2 = nodes[a], nodes[b]
        h = x2 - x1
        # du_h = sum_a u_a * dphi_a/dx, with dphi_a/dx = -1/h, dphi_b/dx = +1/h
        du_h = u_h[a] * (-1.0 / h) + u_h[b] * (1.0 / h)
        total += quadrature(lambda x, du_h=du_h: (du_h - du_ex(x)) ** 2, x1, x2)
    return np.sqrt(total)

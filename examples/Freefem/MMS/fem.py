"""FEM toolbox: quadrature, assembly, error norms.

Two parallel families of routines:

* **1D scalar-integrand** path used by the 1D bar driver:
  `line_quadrature`, `assemble_nodal_forces_1d`, `l2_error_1d`,
  `h1_semi_error_1d`.

* **Dim-agnostic vector path** used by the 2D and 3D drivers:
  `quad_q1_rule`, `tri_p1_rule`, `hex_q1_rule` (element rules),
  `edge_line_rule`, `quad_face_rule` (facet rules),
  `assemble_nodal_forces`, `assemble_traction`,
  `l2_error`, `h1_semi_error`.

Element-rule protocol:

    rule(xe) yields (coords, w, N, dN_phys)
        xe       : (n_local, dim) physical node coordinates of one element
        coords   : (dim,)         physical Gauss-point coordinates
        w        : scalar         weight (already × detJ / triangle area)
        N        : (n_local,)     shape values at the Gauss point
        dN_phys  : (dim, n_local) physical-coord shape gradients, indexed
                                  by (spatial axis, local node)

Facet-rule protocol:

    rule(xe) yields (coords, w, N)
        same convention; no shape gradient (not needed on boundary).
"""
import numpy as np

# ---------------------------------------------------------------------------
# Canonical 1D Gauss-Legendre table on [-1, 1]
# ---------------------------------------------------------------------------

_GAUSS_LEGENDRE_1D = {
    1: (np.array([0.0]),
        np.array([2.0])),
    2: (np.array([-1.0 / np.sqrt(3.0), 1.0 / np.sqrt(3.0)]),
        np.array([1.0, 1.0])),
}


# ---------------------------------------------------------------------------
# 1D scalar-integrand quadrature  (consumed by the 1D bar driver only)
# ---------------------------------------------------------------------------

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
L2_QUADRATURE_1D = line_quadrature(2)
H1_QUADRATURE_1D = line_quadrature(2)


# ---------------------------------------------------------------------------
# Reference-element shape functions
# ---------------------------------------------------------------------------

def _shape_q1(xi, eta):
    """2D Q1 shape on [-1, 1]^2. Node order: (-,-), (+,-), (+,+), (-,+)."""
    N = 0.25 * np.array([
        (1 - xi) * (1 - eta),
        (1 + xi) * (1 - eta),
        (1 + xi) * (1 + eta),
        (1 - xi) * (1 + eta),
    ])
    dN_dxi  = 0.25 * np.array([-(1 - eta),  (1 - eta), (1 + eta), -(1 + eta)])
    dN_deta = 0.25 * np.array([-(1 - xi),  -(1 + xi),  (1 + xi),  (1 - xi) ])
    return N, dN_dxi, dN_deta


def _shape_hex_q1(xi, eta, zeta):
    """3D Q1 hex shape on [-1, 1]^3. Node order matches SOFA RegularGridTopology:
       0:(-,-,-) 1:(+,-,-) 2:(+,+,-) 3:(-,+,-)
       4:(-,-,+) 5:(+,-,+) 6:(+,+,+) 7:(-,+,+)."""
    N = 0.125 * np.array([
        (1 - xi) * (1 - eta) * (1 - zeta),
        (1 + xi) * (1 - eta) * (1 - zeta),
        (1 + xi) * (1 + eta) * (1 - zeta),
        (1 - xi) * (1 + eta) * (1 - zeta),
        (1 - xi) * (1 - eta) * (1 + zeta),
        (1 + xi) * (1 - eta) * (1 + zeta),
        (1 + xi) * (1 + eta) * (1 + zeta),
        (1 - xi) * (1 + eta) * (1 + zeta),
    ])
    dN_dxi = 0.125 * np.array([
        -(1 - eta) * (1 - zeta),  (1 - eta) * (1 - zeta),
         (1 + eta) * (1 - zeta), -(1 + eta) * (1 - zeta),
        -(1 - eta) * (1 + zeta),  (1 - eta) * (1 + zeta),
         (1 + eta) * (1 + zeta), -(1 + eta) * (1 + zeta),
    ])
    dN_deta = 0.125 * np.array([
        -(1 - xi) * (1 - zeta), -(1 + xi) * (1 - zeta),
         (1 + xi) * (1 - zeta),  (1 - xi) * (1 - zeta),
        -(1 - xi) * (1 + zeta), -(1 + xi) * (1 + zeta),
         (1 + xi) * (1 + zeta),  (1 - xi) * (1 + zeta),
    ])
    dN_dzeta = 0.125 * np.array([
        -(1 - xi) * (1 - eta), -(1 + xi) * (1 - eta),
        -(1 + xi) * (1 + eta), -(1 - xi) * (1 + eta),
         (1 - xi) * (1 - eta),  (1 + xi) * (1 - eta),
         (1 + xi) * (1 + eta),  (1 - xi) * (1 + eta),
    ])
    return N, dN_dxi, dN_deta, dN_dzeta


# Reference-triangle quadrature, integrating over {(xi,eta): xi,eta>=0, xi+eta<=1}.
# Weights sum to the reference area 1/2; the physical weight is w_ref * 2*area.
_TRI_QUADRATURE = {
    1: (np.array([[1/3, 1/3]]),
        np.array([1/2])),
    3: (np.array([[1/6, 1/6], [2/3, 1/6], [1/6, 2/3]]),
        np.array([1/6, 1/6, 1/6])),
}


# ---------------------------------------------------------------------------
# Element rules  (consume by assemble_nodal_forces, l2_error, h1_semi_error)
# ---------------------------------------------------------------------------

def quad_q1_rule(n_pts=2):
    """Element rule for Q1 quads: tensor-product Gauss-Legendre on [-1,1]^2."""
    if n_pts not in _GAUSS_LEGENDRE_1D:
        raise ValueError(f"quad_q1_rule: {n_pts}-point rule not supported")
    xi_pts, w_pts = _GAUSS_LEGENDRE_1D[n_pts]

    def rule(xe):
        xe_arr = np.asarray(xe, float)            # (4, 2)
        x_col, y_col = xe_arr[:, 0], xe_arr[:, 1]
        for xi, wi in zip(xi_pts, w_pts):
            for eta, wj in zip(xi_pts, w_pts):
                N, dN_dxi, dN_deta = _shape_q1(xi, eta)
                J = np.array([[dN_dxi  @ x_col, dN_dxi  @ y_col],
                              [dN_deta @ x_col, dN_deta @ y_col]])
                detJ   = np.linalg.det(J)
                Jinv   = np.linalg.inv(J)
                coords = np.array([N @ x_col, N @ y_col])
                dN_dx  = Jinv[0, 0] * dN_dxi + Jinv[1, 0] * dN_deta
                dN_dy  = Jinv[0, 1] * dN_dxi + Jinv[1, 1] * dN_deta
                dN_phys = np.stack([dN_dx, dN_dy])      # (2, 4) [d, a]
                yield coords, wi * wj * detJ, N, dN_phys
    return rule


def tri_p1_rule(n_pts=1):
    """Element rule for P1 triangles: 1-point (centroid) or 3-point Gauss.

    Shape gradients are constant per triangle. Node order in `xe` is the
    canonical P1 ordering — corresponds to (N0, N1, N2) = (1-xi-eta, xi, eta).
    """
    if n_pts not in _TRI_QUADRATURE:
        raise ValueError(f"tri_p1_rule: {n_pts}-point rule not supported")
    pts, wts = _TRI_QUADRATURE[n_pts]

    def rule(xe):
        xe_arr = np.asarray(xe, float)            # (3, 2)
        (x0, y0), (x1, y1), (x2, y2) = xe_arr
        A2    = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
        area  = abs(A2) / 2.0
        dN_dx = np.array([(y1 - y2) / A2, (y2 - y0) / A2, (y0 - y1) / A2])
        dN_dy = np.array([(x2 - x1) / A2, (x0 - x2) / A2, (x1 - x0) / A2])
        dN_phys = np.stack([dN_dx, dN_dy])        # (2, 3) [d, a]
        for (xi, eta), w_ref in zip(pts, wts):
            N      = np.array([1.0 - xi - eta, xi, eta])
            coords = np.array([N @ xe_arr[:, 0], N @ xe_arr[:, 1]])
            yield coords, w_ref * 2.0 * area, N, dN_phys
    return rule


def hex_q1_rule(n_pts=2):
    """Element rule for Q1 hexes: tensor-product Gauss-Legendre on [-1,1]^3."""
    if n_pts not in _GAUSS_LEGENDRE_1D:
        raise ValueError(f"hex_q1_rule: {n_pts}-point rule not supported")
    xi_pts, w_pts = _GAUSS_LEGENDRE_1D[n_pts]

    def rule(xe):
        xe_arr = np.asarray(xe, float)            # (8, 3)
        x_col, y_col, z_col = xe_arr[:, 0], xe_arr[:, 1], xe_arr[:, 2]
        for xi, wi in zip(xi_pts, w_pts):
            for eta, wj in zip(xi_pts, w_pts):
                for zeta, wk in zip(xi_pts, w_pts):
                    N, dN_dxi, dN_deta, dN_dzeta = _shape_hex_q1(xi, eta, zeta)
                    J = np.array([
                        [dN_dxi   @ x_col, dN_dxi   @ y_col, dN_dxi   @ z_col],
                        [dN_deta  @ x_col, dN_deta  @ y_col, dN_deta  @ z_col],
                        [dN_dzeta @ x_col, dN_dzeta @ y_col, dN_dzeta @ z_col],
                    ])
                    detJ   = np.linalg.det(J)
                    Jinv   = np.linalg.inv(J)
                    coords = np.array([N @ x_col, N @ y_col, N @ z_col])
                    dN_dx  = Jinv[0, 0] * dN_dxi + Jinv[1, 0] * dN_deta + Jinv[2, 0] * dN_dzeta
                    dN_dy  = Jinv[0, 1] * dN_dxi + Jinv[1, 1] * dN_deta + Jinv[2, 1] * dN_dzeta
                    dN_dz  = Jinv[0, 2] * dN_dxi + Jinv[1, 2] * dN_deta + Jinv[2, 2] * dN_dzeta
                    dN_phys = np.stack([dN_dx, dN_dy, dN_dz])   # (3, 8) [d, a]
                    yield coords, wi * wj * wk * detJ, N, dN_phys
    return rule


# ---------------------------------------------------------------------------
# Facet rules  (consumed by assemble_traction)
# ---------------------------------------------------------------------------

def edge_line_rule(n_pts=2):
    """Facet rule for a 2D boundary edge: linear shape, 1D Gauss-Legendre.

    Yields per Gauss point along the edge: (coords (2,), w, N (2,)).
    """
    if n_pts not in _GAUSS_LEGENDRE_1D:
        raise ValueError(f"edge_line_rule: {n_pts}-point rule not supported")
    xi_pts, w_pts = _GAUSS_LEGENDRE_1D[n_pts]

    def rule(xe):
        xe_arr = np.asarray(xe, float)            # (2, 2)
        Le = np.linalg.norm(xe_arr[1] - xe_arr[0])
        for xi, wi in zip(xi_pts, w_pts):
            t = 0.5 * (xi + 1.0)
            N = np.array([1.0 - t, t])
            coords = N @ xe_arr                   # (2,)
            yield coords, wi * Le / 2.0, N
    return rule


def quad_face_rule(n_pts=2):
    """Facet rule for a 3D boundary quad face: bilinear shape, 2D Gauss.

    Yields per Gauss point on the face: (coords (3,), w, N (4,)).
    Quad orientation does not affect the magnitude `dA = ||t_xi × t_eta||`.
    """
    if n_pts not in _GAUSS_LEGENDRE_1D:
        raise ValueError(f"quad_face_rule: {n_pts}-point rule not supported")
    xi_pts, w_pts = _GAUSS_LEGENDRE_1D[n_pts]

    def rule(xe):
        xe_arr = np.asarray(xe, float)            # (4, 3)
        for xi, wi in zip(xi_pts, w_pts):
            for eta, wj in zip(xi_pts, w_pts):
                N = 0.25 * np.array([
                    (1 - xi) * (1 - eta),
                    (1 + xi) * (1 - eta),
                    (1 + xi) * (1 + eta),
                    (1 - xi) * (1 + eta),
                ])
                dN_dxi  = 0.25 * np.array([-(1 - eta),  (1 - eta), (1 + eta), -(1 + eta)])
                dN_deta = 0.25 * np.array([-(1 - xi),  -(1 + xi),  (1 + xi),  (1 - xi) ])
                t_xi   = dN_dxi  @ xe_arr         # (3,)
                t_eta  = dN_deta @ xe_arr
                dA     = np.linalg.norm(np.cross(t_xi, t_eta))
                coords = N @ xe_arr               # (3,)
                yield coords, wi * wj * dA, N
    return rule


# ---------------------------------------------------------------------------
# 1D FEM assembly  (scalar-integrand path)
# ---------------------------------------------------------------------------

def assemble_nodal_forces_1d(f_body, nodes, edges, quadrature):
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
# Dim-agnostic FEM assembly  (2D and 3D)
# ---------------------------------------------------------------------------

def assemble_nodal_forces(f_body, nodes, conn, element_rule):
    """
    Assemble consistent nodal body forces, dim- and element-type agnostic.

        F_a = sum_e integral_{Omega_e} N_a(x) * f_body(x) dx

    f_body       : callable (*coords) -> dim-tuple of force components
    nodes        : (N, dim) array
    conn         : iterable of node-index lists, one per element
    element_rule : rule from quad_q1_rule / tri_p1_rule / hex_q1_rule
    """
    nodes = np.asarray(nodes)
    dim   = nodes.shape[1]
    F     = np.zeros((len(nodes), dim))
    for elem in conn:
        # np.asarray forces row-style fancy indexing even when `elem` is a
        # Python tuple (which would otherwise trigger multi-axis indexing).
        idx = np.asarray(elem)
        xe  = nodes[idx]                          # (n_local, dim)
        for coords, w, N, _ in element_rule(xe):
            f_components = f_body(*coords)
            for a, node in enumerate(idx):
                for d in range(dim):
                    F[node, d] += N[a] * f_components[d] * w
    return F


def assemble_traction(traction, nodes, facets, facet_rule):
    """
    Assemble consistent nodal forces from a boundary traction over facets.

        F_a += integral_{facet} N_a(s) * traction(x(s)) dS

    The outward unit normal is baked into `traction` by the caller (one
    lambda per boundary face); the facet rule only supplies the surface area
    element and shape values.

    traction   : callable (*coords) -> dim-tuple of force components
    nodes      : (N, dim) array
    facets     : iterable of node-index tuples (one per boundary facet)
    facet_rule : rule from edge_line_rule (2D) / quad_face_rule (3D)
    """
    nodes = np.asarray(nodes)
    dim   = nodes.shape[1]
    F     = np.zeros((len(nodes), dim))
    for facet in facets:
        idx = np.asarray(facet)
        xe  = nodes[idx]                          # (n_local, dim)
        for coords, w, N in facet_rule(xe):
            T = traction(*coords)
            for a, node in enumerate(idx):
                for d in range(dim):
                    F[node, d] += N[a] * T[d] * w
    return F


# ---------------------------------------------------------------------------
# 1D error norms
# ---------------------------------------------------------------------------

def l2_error_1d(nodes, edges, u_h, u_ex, quadrature):
    """L2 error norm: sqrt( integral (u_h - u_ex)^2 dx ) over the mesh."""
    total = 0.0
    for a, b in edges:
        x1, x2 = nodes[a], nodes[b]
        h = x2 - x1
        u_a, u_b = u_h[a], u_h[b]
        u_interp = lambda x, x1=x1, h=h, u_a=u_a, u_b=u_b: u_a + (u_b - u_a) * (x - x1) / h
        total += quadrature(lambda x: (u_interp(x) - u_ex(x)) ** 2, x1, x2)
    return np.sqrt(total)


def h1_semi_error_1d(nodes, edges, u_h, du_ex, quadrature):
    """H1 semi-norm error: sqrt( integral (du_h - du_ex)^2 dx ) over the mesh."""
    total = 0.0
    for a, b in edges:
        x1, x2 = nodes[a], nodes[b]
        h = x2 - x1
        # du_h = sum_a u_a * dphi_a/dx, with dphi_a/dx = -1/h, dphi_b/dx = +1/h
        du_h = u_h[a] * (-1.0 / h) + u_h[b] * (1.0 / h)
        total += quadrature(lambda x, du_h=du_h: (du_h - du_ex(x)) ** 2, x1, x2)
    return np.sqrt(total)


# ---------------------------------------------------------------------------
# Dim-agnostic error norms (2D and 3D)
# ---------------------------------------------------------------------------

def l2_error(nodes, conn, u_h, u_ex, element_rule):
    """
    Vector L2 error norm, dim- and element-type agnostic.

        sqrt( integral || u_h(x) - u_ex(x) ||^2 dx ) over the mesh

    nodes        : (N, dim) array
    conn         : iterable of node-index lists, one per element
    u_h          : (N, dim) array of nodal displacements
    u_ex         : callable (*coords) -> dim-tuple
    element_rule : rule from quad_q1_rule / tri_p1_rule / hex_q1_rule
    """
    nodes = np.asarray(nodes)
    u     = np.asarray(u_h)
    total = 0.0
    for elem in conn:
        idx   = np.asarray(elem)
        xe    = nodes[idx]                        # (n_local, dim)
        u_loc = u[idx]                            # (n_local, dim)
        for coords, w, N, _ in element_rule(xe):
            u_h_g  = N @ u_loc                    # (dim,)
            u_ex_g = np.asarray(u_ex(*coords))    # (dim,)
            diff   = u_h_g - u_ex_g
            total += float(np.dot(diff, diff)) * w
    return np.sqrt(total)


def h1_semi_error(nodes, conn, u_h, grad_u_ex, element_rule):
    """
    Vector H1 semi-norm error, dim- and element-type agnostic.

        sqrt( integral || grad u_h - grad u_ex ||_F^2 dx ) over the mesh

    nodes        : (N, dim) array
    conn         : iterable of node-index lists, one per element
    u_h          : (N, dim) array of nodal displacements
    grad_u_ex    : callable (*coords) -> (dim, dim) array with [i, j] = du_i/dx_j
    element_rule : rule from quad_q1_rule / tri_p1_rule / hex_q1_rule
    """
    nodes = np.asarray(nodes)
    u     = np.asarray(u_h)
    total = 0.0
    for elem in conn:
        idx   = np.asarray(elem)
        xe    = nodes[idx]                        # (n_local, dim)
        u_loc = u[idx]                            # (n_local, dim)
        for coords, w, _, dN_phys in element_rule(xe):
            # grad_uh[i, d] = sum_a u_i[a] * dN_a/dx_d
            grad_uh = (dN_phys @ u_loc).T         # (dim, dim) [i, d]
            grad_ue = np.asarray(grad_u_ex(*coords))  # [i, j] = du_i/dx_j
            diff    = grad_uh - grad_ue
            total  += float(np.sum(diff * diff)) * w
    return np.sqrt(total)

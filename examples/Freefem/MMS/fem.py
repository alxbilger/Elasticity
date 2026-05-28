"""FEM toolbox: quadrature, assembly, error norms.
"""
import numpy as np

# ---------------------------------------------------------------------------
# Quadrature rules
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
L2_QUADRATURE_1D = line_quadrature(2)
H1_QUADRATURE_1D = line_quadrature(2)


# ---------------------------------------------------------------------------
# 2D Q1 reference shape functions on [-1,1]^2.
# Node ordering: (-1,-1), (+1,-1), (+1,+1), (-1,+1)
# ---------------------------------------------------------------------------

def _shape_q1(xi, eta):
    N = 0.25 * np.array([
        (1 - xi) * (1 - eta),
        (1 + xi) * (1 - eta),
        (1 + xi) * (1 + eta),
        (1 - xi) * (1 + eta),
    ])
    dN_dxi  = 0.25 * np.array([-(1 - eta),  (1 - eta), (1 + eta), -(1 + eta)])
    dN_deta = 0.25 * np.array([-(1 - xi),  -(1 + xi),  (1 + xi),  (1 - xi)])
    return N, dN_dxi, dN_deta


# ---------------------------------------------------------------------------
# 2D element rules
#
# An "element rule" is a callable rule(xe, ye) that yields one tuple
#     (xg, yg, w, N, dN_dx, dN_dy)
# per Gauss point of a single physical element, with:
#   (xg, yg)        physical Gauss point
#   w               weight including detJ (or triangle area)
#   N               length-n_local shape values at the Gauss point
#   dN_dx, dN_dy    length-n_local physical-coord shape gradients
#
# The same protocol is consumed by assemble_nodal_forces_2d, l2_error_2d,
# and h1_semi_error_2d, so adding a new element type is one rule function.
# ---------------------------------------------------------------------------

def quad_q1_rule(n_pts=2):
    """Element rule for Q1 quads: tensor-product Gauss-Legendre."""
    if n_pts not in _GAUSS_LEGENDRE_1D:
        raise ValueError(f"quad_q1_rule: {n_pts}-point rule not supported")
    xi_pts, w_pts = _GAUSS_LEGENDRE_1D[n_pts]

    def rule(xe, ye):
        for xi, wi in zip(xi_pts, w_pts):
            for eta, wj in zip(xi_pts, w_pts):
                N, dN_dxi, dN_deta = _shape_q1(xi, eta)
                J = np.array([[dN_dxi  @ xe, dN_dxi  @ ye],
                              [dN_deta @ xe, dN_deta @ ye]])
                detJ = np.linalg.det(J)
                Jinv = np.linalg.inv(J)
                xg, yg = N @ xe, N @ ye
                dN_dx = Jinv[0, 0] * dN_dxi + Jinv[1, 0] * dN_deta
                dN_dy = Jinv[0, 1] * dN_dxi + Jinv[1, 1] * dN_deta
                yield xg, yg, wi * wj * detJ, N, dN_dx, dN_dy
    return rule


# Reference-triangle quadrature, integrating over {(xi,eta): xi,eta>=0, xi+eta<=1}.
# Weights sum to the reference area 1/2; the physical weight is w_ref * 2*area.
_TRI_QUADRATURE = {
    1: (np.array([[1/3, 1/3]]),
        np.array([1/2])),
    3: (np.array([[1/6, 1/6], [2/3, 1/6], [1/6, 2/3]]),
        np.array([1/6, 1/6, 1/6])),
}


def tri_p1_rule(n_pts=1):
    """Element rule for P1 triangles: 1-point (centroid) or 3-point Gauss.

    Shape gradients are constant per triangle. Node order in (xe, ye) is the
    canonical P1 ordering — corresponds to (N0, N1, N2) = (1-xi-eta, xi, eta).
    """
    if n_pts not in _TRI_QUADRATURE:
        raise ValueError(f"tri_p1_rule: {n_pts}-point rule not supported")
    pts, wts = _TRI_QUADRATURE[n_pts]

    def rule(xe, ye):
        x0, x1, x2 = xe
        y0, y1, y2 = ye
        # Signed double-area: positive for CCW orientation
        A2    = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
        area  = abs(A2) / 2.0
        # Constant physical-coord gradients of the P1 shape functions
        dN_dx = np.array([(y1 - y2) / A2, (y2 - y0) / A2, (y0 - y1) / A2])
        dN_dy = np.array([(x2 - x1) / A2, (x0 - x2) / A2, (x1 - x0) / A2])
        for (xi, eta), w_ref in zip(pts, wts):
            N    = np.array([1.0 - xi - eta, xi, eta])
            xg   = N[0] * x0 + N[1] * x1 + N[2] * x2
            yg   = N[0] * y0 + N[1] * y1 + N[2] * y2
            # w_ref integrates over reference area 1/2; scale by 2*area for physical
            yield xg, yg, w_ref * 2.0 * area, N, dN_dx, dN_dy
    return rule


# ---------------------------------------------------------------------------
# 3D Q1 reference shape functions on [-1,1]^3.
# Node ordering matches SOFA's RegularGridTopology hex convention:
#   0:(-,-,-) 1:(+,-,-) 2:(+,+,-) 3:(-,+,-)
#   4:(-,-,+) 5:(+,-,+) 6:(+,+,+) 7:(-,+,+)
# ---------------------------------------------------------------------------

def _shape_hex_q1(xi, eta, zeta):
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


# ---------------------------------------------------------------------------
# 3D element rules
#
# Same protocol as the 2D rules, lifted one dimension:
#   rule(xe, ye, ze) yields (xg, yg, zg, w, N, dN_dx, dN_dy, dN_dz)
# ---------------------------------------------------------------------------

def hex_q1_rule(n_pts=2):
    """Element rule for Q1 hexes: tensor-product Gauss-Legendre over [-1,1]^3."""
    if n_pts not in _GAUSS_LEGENDRE_1D:
        raise ValueError(f"hex_q1_rule: {n_pts}-point rule not supported")
    xi_pts, w_pts = _GAUSS_LEGENDRE_1D[n_pts]

    def rule(xe, ye, ze):
        for xi, wi in zip(xi_pts, w_pts):
            for eta, wj in zip(xi_pts, w_pts):
                for zeta, wk in zip(xi_pts, w_pts):
                    N, dN_dxi, dN_deta, dN_dzeta = _shape_hex_q1(xi, eta, zeta)
                    J = np.array([
                        [dN_dxi   @ xe, dN_dxi   @ ye, dN_dxi   @ ze],
                        [dN_deta  @ xe, dN_deta  @ ye, dN_deta  @ ze],
                        [dN_dzeta @ xe, dN_dzeta @ ye, dN_dzeta @ ze],
                    ])
                    detJ = np.linalg.det(J)
                    Jinv = np.linalg.inv(J)
                    xg, yg, zg = N @ xe, N @ ye, N @ ze
                    dN_dx = Jinv[0, 0] * dN_dxi + Jinv[1, 0] * dN_deta + Jinv[2, 0] * dN_dzeta
                    dN_dy = Jinv[0, 1] * dN_dxi + Jinv[1, 1] * dN_deta + Jinv[2, 1] * dN_dzeta
                    dN_dz = Jinv[0, 2] * dN_dxi + Jinv[1, 2] * dN_deta + Jinv[2, 2] * dN_dzeta
                    yield xg, yg, zg, wi * wj * wk * detJ, N, dN_dx, dN_dy, dN_dz
    return rule


# ---------------------------------------------------------------------------
# FEM assembly
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


def assemble_traction_2d(traction, nodes, edges, n_pts=2):
    """
    Assemble consistent nodal forces from a boundary traction along edges:

        F_a += integral_{edge} N_a(s) * traction(x(s), y(s)) ds

    Linear edge shape functions in physical arc-length. The outward normal is
    baked into `traction` by the caller (one lambda per boundary side).

    traction : callable (x, y) -> (Tx, Ty)
    nodes    : (N, 2) or (N, 3) array — only the first two columns used
    edges    : iterable of (a, b) node-index pairs along the boundary
    n_pts    : Gauss points per edge (default 2)
    """
    if n_pts not in _GAUSS_LEGENDRE_1D:
        raise ValueError(f"assemble_traction_2d: {n_pts}-point rule not supported")
    xi_pts, w_pts = _GAUSS_LEGENDRE_1D[n_pts]

    xy = np.asarray(nodes)[:, :2]
    F  = np.zeros((len(xy), 2))
    for a, b in edges:
        x1, y1 = xy[a]
        x2, y2 = xy[b]
        Le = np.hypot(x2 - x1, y2 - y1)
        for xi, wi in zip(xi_pts, w_pts):
            t = 0.5 * (xi + 1.0)
            xg = (1.0 - t) * x1 + t * x2
            yg = (1.0 - t) * y1 + t * y2
            Tx, Ty = traction(xg, yg)
            N0, N1 = 1.0 - t, t
            w = wi * Le / 2.0
            F[a, 0] += N0 * Tx * w
            F[a, 1] += N0 * Ty * w
            F[b, 0] += N1 * Tx * w
            F[b, 1] += N1 * Ty * w
    return F


def assemble_nodal_forces_2d(f_body, nodes, conn, element_rule):
    """
    Assemble consistent nodal forces on a 2D mesh, agnostic to element type.

    F_a = sum_e integral_{Omega_e} N_a(x,y) * f_body(x,y) dx dy

    f_body       : callable (x, y) -> (fx, fy)
    nodes        : (N, 2) or (N, 3) array — only the first two columns used
    conn         : iterable of node-index lists, one per element
    element_rule : rule from quad_q1_rule(n_pts) / tri_p1_rule(n_pts)
    """
    xy = np.asarray(nodes)[:, :2]
    F  = np.zeros((len(xy), 2))
    for elem in conn:
        xe, ye = xy[elem, 0], xy[elem, 1]
        for xg, yg, w, N, _, _ in element_rule(xe, ye):
            fx, fy = f_body(xg, yg)
            for a, node in enumerate(elem):
                F[node, 0] += N[a] * fx * w
                F[node, 1] += N[a] * fy * w
    return F


def assemble_traction_3d(traction, nodes, quads, n_pts=2):
    """
    Assemble consistent nodal forces from a boundary traction along quad faces:

        F_a += integral_{face} N_a(s,t) * traction(x(s,t), y(s,t), z(s,t)) dA

    Bilinear quad shape functions on each face. The outward unit normal is
    baked into `traction` by the caller (one lambda per boundary face) — same
    pattern as `assemble_traction_2d`.

    traction : callable (x, y, z) -> (Tx, Ty, Tz)
    nodes    : (N, 3) array
    quads    : iterable of (a, b, c, d) node-index 4-tuples on the boundary
    n_pts    : Gauss points per face direction (default 2 → 2×2 = 4 pts)
    """
    if n_pts not in _GAUSS_LEGENDRE_1D:
        raise ValueError(f"assemble_traction_3d: {n_pts}-point rule not supported")
    xi_pts, w_pts = _GAUSS_LEGENDRE_1D[n_pts]

    xyz = np.asarray(nodes)[:, :3]
    F   = np.zeros((len(xyz), 3))
    for quad in quads:
        xe = xyz[quad, 0]
        ye = xyz[quad, 1]
        ze = xyz[quad, 2]
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
                # Surface area element from the cross product of the two tangent vectors
                t_xi  = np.array([dN_dxi  @ xe, dN_dxi  @ ye, dN_dxi  @ ze])
                t_eta = np.array([dN_deta @ xe, dN_deta @ ye, dN_deta @ ze])
                dA = np.linalg.norm(np.cross(t_xi, t_eta))
                xg, yg, zg = N @ xe, N @ ye, N @ ze
                Tx, Ty, Tz = traction(xg, yg, zg)
                w = wi * wj * dA
                for a, node in enumerate(quad):
                    F[node, 0] += N[a] * Tx * w
                    F[node, 1] += N[a] * Ty * w
                    F[node, 2] += N[a] * Tz * w
    return F


def assemble_nodal_forces_3d(f_body, nodes, conn, element_rule):
    """
    Assemble consistent nodal forces on a 3D mesh, agnostic to element type.

    F_a = sum_e integral_{Omega_e} N_a(x,y,z) * f_body(x,y,z) dx dy dz

    f_body       : callable (x, y, z) -> (fx, fy, fz)
    nodes        : (N, 3) array
    conn         : iterable of node-index lists, one per element
    element_rule : rule from hex_q1_rule(n_pts)
    """
    xyz = np.asarray(nodes)[:, :3]
    F   = np.zeros((len(xyz), 3))
    for elem in conn:
        xe, ye, ze = xyz[elem, 0], xyz[elem, 1], xyz[elem, 2]
        for xg, yg, zg, w, N, _, _, _ in element_rule(xe, ye, ze):
            fx, fy, fz = f_body(xg, yg, zg)
            for a, node in enumerate(elem):
                F[node, 0] += N[a] * fx * w
                F[node, 1] += N[a] * fy * w
                F[node, 2] += N[a] * fz * w
    return F


# ---------------------------------------------------------------------------
# Error norms
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


def l2_error_2d(nodes, conn, u_h, u_ex, element_rule):
    """
    Vector L2 error norm on a 2D mesh, element-type agnostic.

        sqrt( integral || u_h(x,y) - u_ex(x,y) ||^2 dx dy ) over the mesh

    nodes        : (N, 2) or (N, 3) array — only the first two columns used
    conn         : iterable of node-index lists, one per element
    u_h          : (N, 2) array of nodal displacements (ux, uy)
    u_ex         : callable (x, y) -> (ux_ex, uy_ex)
    element_rule : rule from quad_q1_rule / tri_p1_rule
    """
    xy = np.asarray(nodes)[:, :2]
    u  = np.asarray(u_h)
    total = 0.0
    for elem in conn:
        xe, ye = xy[elem, 0], xy[elem, 1]
        u_loc  = u[elem]                       # (n_local, 2)
        for xg, yg, w, N, _, _ in element_rule(xe, ye):
            u_h_g = N @ u_loc                  # (2,)
            ux_ex, uy_ex = u_ex(xg, yg)
            ex = u_h_g[0] - ux_ex
            ey = u_h_g[1] - uy_ex
            total += (ex * ex + ey * ey) * w
    return np.sqrt(total)


def h1_semi_error_2d(nodes, conn, u_h, grad_u_ex, element_rule):
    """
    Vector H1 semi-norm error on a 2D mesh, element-type agnostic.

        sqrt( integral || grad u_h - grad u_ex ||_F^2 dx dy ) over the mesh

    nodes        : (N, 2) or (N, 3) array — only the first two columns used
    conn         : iterable of node-index lists, one per element
    u_h          : (N, 2) array of nodal displacements (ux, uy)
    grad_u_ex    : callable (x, y) -> 2x2 array
                   [[dux/dx, dux/dy], [duy/dx, duy/dy]]
    element_rule : rule from quad_q1_rule / tri_p1_rule
    """
    xy = np.asarray(nodes)[:, :2]
    u  = np.asarray(u_h)
    total = 0.0
    for elem in conn:
        xe, ye = xy[elem, 0], xy[elem, 1]
        u_loc  = u[elem]                       # (n_local, 2)
        for xg, yg, w, _, dN_dx, dN_dy in element_rule(xe, ye):
            dux_dx_h = dN_dx @ u_loc[:, 0]
            dux_dy_h = dN_dy @ u_loc[:, 0]
            duy_dx_h = dN_dx @ u_loc[:, 1]
            duy_dy_h = dN_dy @ u_loc[:, 1]
            G = grad_u_ex(xg, yg)
            total += (
                (dux_dx_h - G[0, 0])**2 + (dux_dy_h - G[0, 1])**2 +
                (duy_dx_h - G[1, 0])**2 + (duy_dy_h - G[1, 1])**2
            ) * w
    return np.sqrt(total)


def l2_error_3d(nodes, conn, u_h, u_ex, element_rule):
    """
    Vector L2 error norm on a 3D mesh, element-type agnostic.

        sqrt( integral || u_h(x,y,z) - u_ex(x,y,z) ||^2 dx dy dz ) over the mesh

    nodes        : (N, 3) array
    conn         : iterable of node-index lists, one per element
    u_h          : (N, 3) array of nodal displacements (ux, uy, uz)
    u_ex         : callable (x, y, z) -> (ux_ex, uy_ex, uz_ex)
    element_rule : rule from hex_q1_rule
    """
    xyz = np.asarray(nodes)[:, :3]
    u   = np.asarray(u_h)
    total = 0.0
    for elem in conn:
        xe, ye, ze = xyz[elem, 0], xyz[elem, 1], xyz[elem, 2]
        u_loc = u[elem]                        # (n_local, 3)
        for xg, yg, zg, w, N, _, _, _ in element_rule(xe, ye, ze):
            u_h_g = N @ u_loc                  # (3,)
            ux_ex, uy_ex, uz_ex = u_ex(xg, yg, zg)
            ex = u_h_g[0] - ux_ex
            ey = u_h_g[1] - uy_ex
            ez = u_h_g[2] - uz_ex
            total += (ex * ex + ey * ey + ez * ez) * w
    return np.sqrt(total)


def h1_semi_error_3d(nodes, conn, u_h, grad_u_ex, element_rule):
    """
    Vector H1 semi-norm error on a 3D mesh, element-type agnostic.

        sqrt( integral || grad u_h - grad u_ex ||_F^2 dx dy dz ) over the mesh

    nodes        : (N, 3) array
    conn         : iterable of node-index lists, one per element
    u_h          : (N, 3) array of nodal displacements (ux, uy, uz)
    grad_u_ex    : callable (x, y, z) -> 3x3 array
                   [[dux/dx, dux/dy, dux/dz], [duy/dx, ...], [duz/dx, ...]]
    element_rule : rule from hex_q1_rule
    """
    xyz = np.asarray(nodes)[:, :3]
    u   = np.asarray(u_h)
    total = 0.0
    for elem in conn:
        xe, ye, ze = xyz[elem, 0], xyz[elem, 1], xyz[elem, 2]
        u_loc = u[elem]                        # (n_local, 3)
        for xg, yg, zg, w, _, dN_dx, dN_dy, dN_dz in element_rule(xe, ye, ze):
            dux_dx_h = dN_dx @ u_loc[:, 0]
            dux_dy_h = dN_dy @ u_loc[:, 0]
            dux_dz_h = dN_dz @ u_loc[:, 0]
            duy_dx_h = dN_dx @ u_loc[:, 1]
            duy_dy_h = dN_dy @ u_loc[:, 1]
            duy_dz_h = dN_dz @ u_loc[:, 1]
            duz_dx_h = dN_dx @ u_loc[:, 2]
            duz_dy_h = dN_dy @ u_loc[:, 2]
            duz_dz_h = dN_dz @ u_loc[:, 2]
            G = grad_u_ex(xg, yg, zg)
            total += (
                (dux_dx_h - G[0, 0])**2 + (dux_dy_h - G[0, 1])**2 + (dux_dz_h - G[0, 2])**2 +
                (duy_dx_h - G[1, 0])**2 + (duy_dy_h - G[1, 1])**2 + (duy_dz_h - G[1, 2])**2 +
                (duz_dx_h - G[2, 0])**2 + (duz_dy_h - G[2, 1])**2 + (duz_dz_h - G[2, 2])**2
            ) * w
    return np.sqrt(total)

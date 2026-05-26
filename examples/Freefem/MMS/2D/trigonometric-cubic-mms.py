"""
MMS 2D - elasticity linear 
2 MMS in this code : 
cubic :   ux = x(L-x)(x+y), uy = y(L-y)(y-x)
trigo :   ux = sin(πx/L) cos(πy/L), uy = cos(πx/L) sin(πy/L)

Q1 & P1 elements 
nu = Poisson ratio 
dim="2d"  ==>  Vec2d  ==>  plane stress   \lambda  = E*nu/(1-nu^2)
dim="3d"  ==>  Vec3d  ==>  plane strain   \lambda  = E* nu/((1+nu)(1-2*nu))
"""

import numpy as np
import matplotlib.pyplot as plt
import os

import Sofa
import Sofa.Core
import Sofa.Simulation
import SofaRuntime

#  =========== Choose the MMS ============================
MMS_TYPE = "cubic"   # "cubic" ou "trigonometric"

RESULTS_DIR = f"results_{MMS_TYPE}_q1p1"
os.makedirs(RESULTS_DIR, exist_ok=True)

GAUSS_PTS = np.array([-1.0 / np.sqrt(3), 1.0 / np.sqrt(3)])
GAUSS_WTS = np.array([1.0, 1.0])

# ============== LAME ====================
def lame(E, nu, dim="2d"):
    """
    dim="2d" ===> plane stress  :  \lambda  = E*nu/(1-nu^²)
    dim="3d"  ==>  Vec3d  ==>  plane strain   \lambda  = E* nu/((1+nu)(1-2*nu))

    """
    if dim == "2d":
        lam = E * nu / (1 - nu**2)
    else:
        lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    return lam, mu

# ==================== Cubic MMS =================================
def ux_mms_cubic(x, y, L):
    return x * (L - x) * (x + y)

def uy_mms_cubic(x, y, L):
    return y * (L - y) * (y - x)

# --------------- Derivatives cubic ---------------
def dux_dx_cubic(x, y, L):
    return (L - 2*x) * (x + y) + x * (L - x)

def dux_dy_cubic(x, y, L):
    return x * (L - x)

def duy_dx_cubic(x, y, L):
    return -y * (L - y)

def duy_dy_cubic(x, y, L):
    return (L - 2*y) * (y - x) + y * (L - y)

def d2ux_dxx_cubic(x, y, L):
    return 2*L - 6*x - 2*y

def d2ux_dxy_cubic(x, y, L):
    return L - 2*x

def d2ux_dyy_cubic(x, y, L):
    return np.zeros_like(np.asarray(x, float))

def d2uy_dxx_cubic(x, y, L):
    return np.zeros_like(np.asarray(x, float))

def d2uy_dxy_cubic(x, y, L):
    return -(L - 2*y)

def d2uy_dyy_cubic(x, y, L):
    return 2*L - 6*y + 2*x

# --------------- Body force cubic ---------------
def fx_body_cubic(x, y, E, nu, L, dim="2d"):
    lam, mu = lame(E, nu, dim)
    return -((lam + 2*mu) * d2ux_dxx_cubic(x, y, L)
             + lam        * d2uy_dxy_cubic(x, y, L)
             + mu         * (d2ux_dyy_cubic(x, y, L) + d2uy_dxy_cubic(x, y, L)))

def fy_body_cubic(x, y, E, nu, L, dim="2d"):
    lam, mu = lame(E, nu, dim)
    return -(mu          * (d2ux_dxy_cubic(x, y, L) + d2uy_dxx_cubic(x, y, L))
             + lam        * d2ux_dxy_cubic(x, y, L)
             + (lam + 2*mu) * d2uy_dyy_cubic(x, y, L))

# --------------- Stress cubic ---------------
def sigma_xx_cubic(x, y, E, nu, L, dim="2d"):
    lam, mu = lame(E, nu, dim)
    return (lam + 2*mu) * dux_dx_cubic(x, y, L) + lam * duy_dy_cubic(x, y, L)

def sigma_yy_cubic(x, y, E, nu, L, dim="2d"):
    lam, mu = lame(E, nu, dim)
    return lam * dux_dx_cubic(x, y, L) + (lam + 2*mu) * duy_dy_cubic(x, y, L)

def sigma_xy_cubic(x, y, E, nu, L, dim="2d"):
    _, mu = lame(E, nu, dim)
    return mu * (dux_dy_cubic(x, y, L) + duy_dx_cubic(x, y, L))

def traction_cubic(x, y, nx_c, ny_c, E, nu, L, dim="2d"):
    sxx = sigma_xx_cubic(x, y, E, nu, L, dim)
    syy = sigma_yy_cubic(x, y, E, nu, L, dim)
    sxy = sigma_xy_cubic(x, y, E, nu, L, dim)
    return sxx*nx_c + sxy*ny_c, sxy*nx_c + syy*ny_c

# ==================== Trigonometric MMS ===================================
def ux_mms_trigo(x, y, L):
    return np.sin(np.pi * x / L) * np.cos(np.pi * y / L)

def uy_mms_trigo(x, y, L):
    return np.cos(np.pi * x / L) * np.sin(np.pi * y / L)

# --------------- Trigo derivatives ---------------
def dux_dx_trigo(x, y, L):
    return (np.pi / L) * np.cos(np.pi * x / L) * np.cos(np.pi * y / L)

def dux_dy_trigo(x, y, L):
    return -(np.pi / L) * np.sin(np.pi * x / L) * np.sin(np.pi * y / L)

def duy_dx_trigo(x, y, L):
    return -(np.pi / L) * np.sin(np.pi * x / L) * np.sin(np.pi * y / L)

def duy_dy_trigo(x, y, L):
    return (np.pi / L) * np.cos(np.pi * x / L) * np.cos(np.pi * y / L)

def d2ux_dxx_trigo(x, y, L):
    return -(np.pi / L)**2 * ux_mms_trigo(x, y, L)

def d2ux_dyy_trigo(x, y, L):
    return -(np.pi / L)**2 * ux_mms_trigo(x, y, L)

def d2ux_dxy_trigo(x, y, L):
    return -(np.pi / L)**2 * np.cos(np.pi * x / L) * np.sin(np.pi * y / L)

def d2uy_dxx_trigo(x, y, L):
    return -(np.pi / L)**2 * uy_mms_trigo(x, y, L)

def d2uy_dyy_trigo(x, y, L):
    return -(np.pi / L)**2 * uy_mms_trigo(x, y, L)

def d2uy_dxy_trigo(x, y, L):
    return -(np.pi / L)**2 * np.sin(np.pi * x / L) * np.cos(np.pi * y / L)

# --------------- Trigo body force ---------------
def fx_body_trigo(x, y, E, nu, L, dim="2d"):
    lam, mu = lame(E, nu, dim)
    ddiv_dx  = d2ux_dxx_trigo(x, y, L) + d2uy_dxy_trigo(x, y, L)
    laplace_ux = d2ux_dxx_trigo(x, y, L) + d2ux_dyy_trigo(x, y, L)
    return -(lam + mu) * ddiv_dx - mu * laplace_ux

def fy_body_trigo(x, y, E, nu, L, dim="2d"):
    lam, mu = lame(E, nu, dim)
    ddiv_dy  = d2ux_dxy_trigo(x, y, L) + d2uy_dyy_trigo(x, y, L)
    laplace_uy = d2uy_dxx_trigo(x, y, L) + d2uy_dyy_trigo(x, y, L)
    return -(lam + mu) * ddiv_dy - mu * laplace_uy

# --------------- Stress trigo ---------------
def sigma_xx_trigo(x, y, E, nu, L, dim="2d"):
    lam, mu = lame(E, nu, dim)
    return (lam + 2*mu) * dux_dx_trigo(x, y, L) + lam * duy_dy_trigo(x, y, L)

def sigma_yy_trigo(x, y, E, nu, L, dim="2d"):
    lam, mu = lame(E, nu, dim)
    return lam * dux_dx_trigo(x, y, L) + (lam + 2*mu) * duy_dy_trigo(x, y, L)

def sigma_xy_trigo(x, y, E, nu, L, dim="2d"):
    _, mu = lame(E, nu, dim)
    return mu * (dux_dy_trigo(x, y, L) + duy_dx_trigo(x, y, L))

def traction_trigo(x, y, nx_c, ny_c, E, nu, L, dim="2d"):
    sxx = sigma_xx_trigo(x, y, E, nu, L, dim)
    syy = sigma_yy_trigo(x, y, E, nu, L, dim)
    sxy = sigma_xy_trigo(x, y, E, nu, L, dim)
    return sxx*nx_c + sxy*ny_c, sxy*nx_c + syy*ny_c

# ====================  MMS Selection ====================
if MMS_TYPE == "cubic":
    ux_mms = ux_mms_cubic
    uy_mms = uy_mms_cubic
    dux_dx = dux_dx_cubic
    dux_dy = dux_dy_cubic
    duy_dx = duy_dx_cubic
    duy_dy = duy_dy_cubic
    fx_body  = fx_body_cubic
    fy_body  = fy_body_cubic
    sigma_xx = sigma_xx_cubic
    sigma_yy = sigma_yy_cubic
    sigma_xy = sigma_xy_cubic
    traction = traction_cubic
    print("   cubic MMS activated")
elif MMS_TYPE == "trigonometric":
    ux_mms = ux_mms_trigo
    uy_mms = uy_mms_trigo
    dux_dx = dux_dx_trigo
    dux_dy = dux_dy_trigo
    duy_dx = duy_dx_trigo
    duy_dy = duy_dy_trigo
    fx_body  = fx_body_trigo
    fy_body  = fy_body_trigo
    sigma_xx = sigma_xx_trigo
    sigma_yy = sigma_yy_trigo
    sigma_xy = sigma_xy_trigo
    traction = traction_trigo
    print("   trigonometric MMS activated")
else:
    raise ValueError(f"MMS_TYPE inconnu: {MMS_TYPE}")

# ==================== MESH ====================
def get_nodes_2d(L, nx, ny, dim="2d"):
    dx, dy = L / (nx - 1), L / (ny - 1)
    if dim == "3d":
        pts = [[i*dx, j*dy, 0.0] for j in range(ny) for i in range(nx)]
    else:
        pts = [[i*dx, j*dy] for j in range(ny) for i in range(nx)]
    return np.array(pts, dtype=float)

def _sofa_template(dim):
    return "Vec3d" if dim == "3d" else "Vec2d"

# ==================== DIRICHLET ====================
def _add_dirichlet(Beam, nodes_2d, L, dim, with_visual):
    tmpl = _sofa_template(dim)
    xy   = nodes_2d[:, :2]
    eps  = 1e-10

    idx_ux0 = [k for k, (xk, _) in enumerate(xy) if xk < eps or xk > L - eps]
    idx_uy0 = [k for k, (_, yk) in enumerate(xy) if yk < eps or yk > L - eps]

    if dim == "2d":
        Beam.addObject("PartialFixedProjectiveConstraint",
                       name="fix_ux", template=tmpl,
                       indices=" ".join(map(str, idx_ux0)),
                       fixedDirections="1 0")
        Beam.addObject("PartialFixedProjectiveConstraint",
                       name="fix_uy", template=tmpl,
                       indices=" ".join(map(str, idx_uy0)),
                       fixedDirections="0 1")
    else:
        all_idx = " ".join(str(k) for k in range(len(nodes_2d)))
        Beam.addObject("PartialFixedProjectiveConstraint",
                       name="fix_ux", template=tmpl,
                       indices=" ".join(map(str, idx_ux0)),
                       fixedDirections="1 0 0")
        Beam.addObject("PartialFixedProjectiveConstraint",
                       name="fix_uy", template=tmpl,
                       indices=" ".join(map(str, idx_uy0)),
                       fixedDirections="0 1 0")
        Beam.addObject("PartialFixedProjectiveConstraint",
                       name="fix_z", template=tmpl,
                       indices=all_idx,
                       fixedDirections="0 0 1")

def _pack_forces(F_xy, dim, n_nodes):
    idx = " ".join(str(k) for k in range(n_nodes))
    if dim == "3d":
        frc = " ".join(f"{F_xy[k,0]} {F_xy[k,1]} 0.0" for k in range(n_nodes))
    else:
        frc = " ".join(f"{F_xy[k,0]} {F_xy[k,1]}" for k in range(n_nodes))
    return idx, frc

# ====================  Q1 element ====================
class _QuadElement:
    LABEL = "Q1 quad"

    @staticmethod
    def get_connectivity(nx, ny):
        quads = []
        for j in range(ny - 1):
            for i in range(nx - 1):
                k00 =  j      * nx + i
                k10 =  j      * nx + (i + 1)
                k01 = (j + 1) * nx + i
                k11 = (j + 1) * nx + (i + 1)
                quads.append([k00, k10, k11, k01])
        return np.array(quads)

    @staticmethod
    def to_triangles(conn):
        tris = []
        for q in conn:
            tris.append([q[0], q[1], q[2]])
            tris.append([q[0], q[2], q[3]])
        return np.array(tris)

    @staticmethod
    def _shape(xi, eta):
        N = 0.25 * np.array([(1-xi)*(1-eta), (1+xi)*(1-eta),
                              (1+xi)*(1+eta), (1-xi)*(1+eta)])
        dN_dxi  = 0.25 * np.array([-(1-eta),  (1-eta),  (1+eta), -(1+eta)])
        dN_deta = 0.25 * np.array([-(1-xi),  -(1+xi),   (1+xi),   (1-xi)])
        return N, dN_dxi, dN_deta

    @staticmethod
    def compute_nodal_forces(nodes_2d, L, E, nu, nx, ny, dim="2d"):
        quads = _QuadElement.get_connectivity(nx, ny)
        xy    = nodes_2d[:, :2]
        F     = np.zeros((len(xy), 2))

        
        for quad in quads:
            xe, ye = xy[quad, 0], xy[quad, 1]
            for xi, wi in zip(GAUSS_PTS, GAUSS_WTS):
                for eta, wj in zip(GAUSS_PTS, GAUSS_WTS):
                    N, dN_dxi, dN_deta = _QuadElement._shape(xi, eta)
                    J = np.array([[dN_dxi @ xe,  dN_dxi @ ye],
                                  [dN_deta @ xe, dN_deta @ ye]])
                    detJ = np.linalg.det(J)
                    xg, yg = N @ xe, N @ ye
                    
                    fx = fx_body(xg, yg, E, nu, L, dim)
                    fy = fy_body(xg, yg, E, nu, L, dim)
                    w = wi * wj * detJ
                    for a, node in enumerate(quad):
                        F[node, 0] += N[a] * fx * w
                        F[node, 1] += N[a] * fy * w

        # Neumann x=0
        for j in range(ny - 1):
            k0, k1 = j * nx, (j + 1) * nx
            x0, y0 = xy[k0]; x1, y1 = xy[k1]
            Le = np.hypot(x1-x0, y1-y0)
            for xi, wi in zip(GAUSS_PTS, GAUSS_WTS):
                t = 0.5 * (xi + 1)
                xg, yg = (1-t)*x0 + t*x1, (1-t)*y0 + t*y1
                Tx, Ty = traction(xg, yg, -1., 0., E, nu, L, dim)
                N0, N1 = 0.5*(1 - xi), 0.5*(1 + xi)
                w = wi * Le / 2.0
                F[k0, 0] += N0 * Tx * w;  F[k0, 1] += N0 * Ty * w
                F[k1, 0] += N1 * Tx * w;  F[k1, 1] += N1 * Ty * w

        # Neumann x=L
        for j in range(ny - 1):
            k0 = j * nx + (nx - 1)
            k1 = (j + 1) * nx + (nx - 1)
            x0, y0 = xy[k0]; x1, y1 = xy[k1]
            Le = np.hypot(x1-x0, y1-y0)
            for xi, wi in zip(GAUSS_PTS, GAUSS_WTS):
                t = 0.5 * (xi + 1)
                xg, yg = (1-t)*x0 + t*x1, (1-t)*y0 + t*y1
                Tx, Ty = traction(xg, yg, +1., 0., E, nu, L, dim)
                N0, N1 = 0.5*(1 - xi), 0.5*(1 + xi)
                w = wi * Le / 2.0
                F[k0, 0] += N0 * Tx * w;  F[k0, 1] += N0 * Ty * w
                F[k1, 0] += N1 * Tx * w;  F[k1, 1] += N1 * Ty * w

        # Neumann y=0
        for i in range(nx - 1):
            k0, k1 = i, i + 1
            x0, y0 = xy[k0]; x1, y1 = xy[k1]
            Le = np.hypot(x1-x0, y1-y0)
            for xi, wi in zip(GAUSS_PTS, GAUSS_WTS):
                t = 0.5 * (xi + 1)
                xg, yg = (1-t)*x0 + t*x1, (1-t)*y0 + t*y1
                Tx, Ty = traction(xg, yg, 0., -1., E, nu, L, dim)
                N0, N1 = 0.5*(1 - xi), 0.5*(1 + xi)
                w = wi * Le / 2.0
                F[k0, 0] += N0 * Tx * w;  F[k0, 1] += N0 * Ty * w
                F[k1, 0] += N1 * Tx * w;  F[k1, 1] += N1 * Ty * w

        # Neumann y=L
        for i in range(nx - 1):
            k0 = (ny - 1) * nx + i
            k1 = (ny - 1) * nx + i + 1
            x0, y0 = xy[k0]; x1, y1 = xy[k1]
            Le = np.hypot(x1-x0, y1-y0)
            for xi, wi in zip(GAUSS_PTS, GAUSS_WTS):
                t = 0.5 * (xi + 1)
                xg, yg = (1-t)*x0 + t*x1, (1-t)*y0 + t*y1
                Tx, Ty = traction(xg, yg, 0., +1., E, nu, L, dim)
                N0, N1 = 0.5*(1 - xi), 0.5*(1 + xi)
                w = wi * Le / 2.0
                F[k0, 0] += N0 * Tx * w;  F[k0, 1] += N0 * Ty * w
                F[k1, 0] += N1 * Tx * w;  F[k1, 1] += N1 * Ty * w

        return F

    @staticmethod
    def create_sofa_scene(rootNode, L=1.0, E=1e6, nu=0.3,
                          nx=10, ny=10, with_visual=True, dim="2d"):
        tmpl = _sofa_template(dim)

        rootNode.addObject("RequiredPlugin", name="Sofa.Component.Visual")
        rootNode.addObject("RequiredPlugin", pluginName=[
            "Elasticity", "Sofa.Component.Constraint.Projective",
            "Sofa.Component.Engine.Select", "Sofa.Component.LinearSolver.Direct",
            "Sofa.Component.MechanicalLoad", "Sofa.Component.ODESolver.Backward",
            "Sofa.Component.StateContainer",
            "Sofa.Component.Topology.Container.Dynamic",
        ])
        rootNode.addObject("DefaultAnimationLoop")
        if with_visual:
            rootNode.addObject("VisualStyle",
                               displayFlags="showBehaviorModels showForceFields")

        nodes_2d = get_nodes_2d(L, nx, ny, dim=dim)
        quads_np = _QuadElement.get_connectivity(nx, ny)

        Beam = rootNode.addChild("Beam")
        Beam.addObject("StaticSolver", name="staticSolver", printLog=False)
        Beam.addObject("SparseLDLSolver", name="linearSolver",
                       template="CompressedRowSparseMatrixd")
        dofs = Beam.addObject("MechanicalObject", name="dofs", template=tmpl,
                              position=nodes_2d.tolist(),
                              showObject=with_visual, showObjectScale=0.005*L)

        Beam.addObject("QuadSetTopologyContainer",
                       name="topology", quads=quads_np.tolist())
        Beam.addObject("QuadSetTopologyModifier")
        Beam.addObject("LinearSmallStrainFEMForceField",
                       name="FEM", template=tmpl,
                       youngModulus=E, poissonRatio=nu, topology="@topology")

        _add_dirichlet(Beam, nodes_2d, L, dim, with_visual)

        
        F_xy = _QuadElement.compute_nodal_forces(nodes_2d, L, E, nu, nx, ny, dim)
        idx, frc = _pack_forces(F_xy, dim, len(nodes_2d))
        Beam.addObject("ConstantForceField", name="MMS_forces", template=tmpl,
                       indices=idx, forces=frc)

        return dofs, nodes_2d, quads_np

    @staticmethod
    def compute_l2(nodes_2d, ux, uy, L, nx, ny, conn):
        xy, err2 = nodes_2d[:, :2], 0.0
        for quad in conn:
            xe, ye = xy[quad, 0], xy[quad, 1]
            for xi, wi in zip(GAUSS_PTS, GAUSS_WTS):
                for eta, wj in zip(GAUSS_PTS, GAUSS_WTS):
                    N, dN_dxi, dN_deta = _QuadElement._shape(xi, eta)
                    J = np.array([[dN_dxi @ xe,  dN_dxi @ ye],
                                  [dN_deta @ xe, dN_deta @ ye]])
                    detJ = np.linalg.det(J)
                    xg, yg = N @ xe, N @ ye
                    ux_h = N @ ux[quad]
                    uy_h = N @ uy[quad]
                    err2 += ((ux_h - ux_mms(xg, yg, L))**2
                           + (uy_h - uy_mms(xg, yg, L))**2) * wi * wj * detJ
        return np.sqrt(err2)

    @staticmethod
    def compute_h1(nodes_2d, ux, uy, L, conn, **_):
        xy, err2 = nodes_2d[:, :2], 0.0
        for quad in conn:
            xe, ye = xy[quad, 0], xy[quad, 1]
            for xi, wi in zip(GAUSS_PTS, GAUSS_WTS):
                for eta, wj in zip(GAUSS_PTS, GAUSS_WTS):
                    N, dN_dxi, dN_deta = _QuadElement._shape(xi, eta)
                    J = np.array([[dN_dxi @ xe,  dN_dxi @ ye],
                                  [dN_deta @ xe, dN_deta @ ye]])
                    detJ = np.linalg.det(J)
                    Jinv = np.linalg.inv(J)
                    xg, yg = N @ xe, N @ ye
                    dN_dx = Jinv[0, 0]*dN_dxi + Jinv[1, 0]*dN_deta
                    dN_dy = Jinv[0, 1]*dN_dxi + Jinv[1, 1]*dN_deta
                    err2 += (
                        (dN_dx @ ux[quad] - dux_dx(xg, yg, L))**2 +
                        (dN_dy @ ux[quad] - dux_dy(xg, yg, L))**2 +
                        (dN_dx @ uy[quad] - duy_dx(xg, yg, L))**2 +
                        (dN_dy @ uy[quad] - duy_dy(xg, yg, L))**2
                    ) * wi * wj * detJ
        return np.sqrt(err2)

# ====================  P1 elements ====================
class _TriElement:
    LABEL = "P1 tri"

    @staticmethod
    def get_connectivity(nx, ny):
        tris = []
        for j in range(ny - 1):
            for i in range(nx - 1):
                k00 =  j      * nx + i
                k10 =  j      * nx + (i + 1)
                k01 = (j + 1) * nx + i
                k11 = (j + 1) * nx + (i + 1)
                tris += [[k00, k10, k11], [k00, k11, k01]]
        return np.array(tris)

    @staticmethod
    def to_triangles(conn):
        return conn

    @staticmethod
    def compute_nodal_forces(nodes_2d, L, E, nu, nx, ny, dim="2d"):
        dx, dy = L / (nx - 1), L / (ny - 1)
        xy = nodes_2d[:, :2]
        F = np.zeros((len(xy), 2))
        eps = 1e-10

    
        for k, (xk, yk) in enumerate(xy):
            i = round(xk / dx)
            j = round(yk / dy)
            wx = dx/2 if (i == 0 or i == nx-1) else dx
            wy = dy/2 if (j == 0 or j == ny-1) else dy
            # dim transmis
            F[k, 0] += fx_body(xk, yk, E, nu, L, dim) * wx * wy
            F[k, 1] += fy_body(xk, yk, E, nu, L, dim) * wx * wy

        # Neumann x=0
        for j in range(ny):
            k = j * nx
            xk, yk = xy[k]
            if yk < eps or yk > L - eps:
                continue
            wy = dy/2 if (j == 0 or j == ny-1) else dy
            Tx, Ty = traction(xk, yk, -1., 0., E, nu, L, dim)
            F[k, 0] += Tx * wy;  F[k, 1] += Ty * wy

        # Neumann x=L
        for j in range(ny):
            k = j * nx + (nx - 1)
            xk, yk = xy[k]
            if yk < eps or yk > L - eps:
                continue
            wy = dy/2 if (j == 0 or j == ny-1) else dy
            Tx, Ty = traction(xk, yk, +1., 0., E, nu, L, dim)
            F[k, 0] += Tx * wy;  F[k, 1] += Ty * wy

        # Neumann y=0
        for i in range(nx):
            k = i
            xk, yk = xy[k]
            if xk < eps or xk > L - eps:
                continue
            wx = dx/2 if (i == 0 or i == nx-1) else dx
            Tx, Ty = traction(xk, yk, 0., -1., E, nu, L, dim)
            F[k, 0] += Tx * wx;  F[k, 1] += Ty * wx

        # Neumann y=L
        for i in range(nx):
            k = (ny - 1) * nx + i
            xk, yk = xy[k]
            if xk < eps or xk > L - eps:
                continue
            wx = dx/2 if (i == 0 or i == nx-1) else dx
            Tx, Ty = traction(xk, yk, 0., +1., E, nu, L, dim)
            F[k, 0] += Tx * wx;  F[k, 1] += Ty * wx

        return F

    @staticmethod
    def create_sofa_scene(rootNode, L=1.0, E=1e6, nu=0.3,
                          nx=10, ny=10, with_visual=True, dim="2d"):
        tmpl = _sofa_template(dim)

        rootNode.addObject("RequiredPlugin", name="Sofa.Component.Visual")
        rootNode.addObject("RequiredPlugin", pluginName=[
            "Elasticity", "Sofa.Component.Constraint.Projective",
            "Sofa.Component.Engine.Select", "Sofa.Component.LinearSolver.Direct",
            "Sofa.Component.MechanicalLoad", "Sofa.Component.ODESolver.Backward",
            "Sofa.Component.StateContainer",
            "Sofa.Component.Topology.Container.Dynamic",
        ])
        rootNode.addObject("DefaultAnimationLoop")
        if with_visual:
            rootNode.addObject("VisualStyle",
                               displayFlags="showBehaviorModels showForceFields")

        nodes_2d = get_nodes_2d(L, nx, ny, dim=dim)
        tris_np = _TriElement.get_connectivity(nx, ny)

        Beam = rootNode.addChild("Beam")
        Beam.addObject("StaticSolver", name="staticSolver", printLog=False)
        Beam.addObject("SparseLDLSolver", name="linearSolver",
                       template="CompressedRowSparseMatrixd")
        dofs = Beam.addObject("MechanicalObject", name="dofs", template=tmpl,
                              position=nodes_2d.tolist(),
                              showObject=with_visual, showObjectScale=0.005*L)

        Beam.addObject("TriangleSetTopologyContainer",
                       name="topology", triangles=tris_np.tolist())
        Beam.addObject("TriangleSetTopologyModifier")
        Beam.addObject("LinearSmallStrainFEMForceField",
                       name="FEM", template=tmpl,
                       youngModulus=E, poissonRatio=nu, topology="@topology")

        _add_dirichlet(Beam, nodes_2d, L, dim, with_visual)

        
        F_xy = _TriElement.compute_nodal_forces(nodes_2d, L, E, nu, nx, ny, dim)
        idx, frc = _pack_forces(F_xy, dim, len(nodes_2d))
        Beam.addObject("ConstantForceField", name="MMS_forces", template=tmpl,
                       indices=idx, forces=frc)

        return dofs, nodes_2d, tris_np

    @staticmethod
    def compute_l2(nodes_2d, ux, uy, L, nx, ny, conn):
        dx, dy = L / (nx - 1), L / (ny - 1)
        area = dx * dy / 2.0
        xy = nodes_2d[:, :2]
        err2 = 0.0
        for i0, i1, i2 in conn:
            xc = (xy[i0,0] + xy[i1,0] + xy[i2,0]) / 3
            yc = (xy[i0,1] + xy[i1,1] + xy[i2,1]) / 3
            ux_c = (ux[i0] + ux[i1] + ux[i2]) / 3
            uy_c = (uy[i0] + uy[i1] + uy[i2]) / 3
            err2 += ((ux_c - ux_mms(xc, yc, L))**2
                   + (uy_c - uy_mms(xc, yc, L))**2) * area
        return np.sqrt(err2)

    @staticmethod
    def _grad_p1(xy, u, tri):
        i0, i1, i2 = tri
        x0, y0 = xy[i0]
        x1, y1 = xy[i1]
        x2, y2 = xy[i2]
        A2 = (x1-x0)*(y2-y0) - (x2-x0)*(y1-y0)
        dudx = (u[i0]*(y1-y2) + u[i1]*(y2-y0) + u[i2]*(y0-y1)) / A2
        dudy = (u[i0]*(x2-x1) + u[i1]*(x0-x2) + u[i2]*(x1-x0)) / A2
        return dudx, dudy, abs(A2) / 2.0

    @staticmethod
    def compute_h1(nodes_2d, ux, uy, L, conn, **_):
        xy, err2 = nodes_2d[:, :2], 0.0
        for tri in conn:
            i0, i1, i2 = tri
            xc = (xy[i0,0] + xy[i1,0] + xy[i2,0]) / 3
            yc = (xy[i0,1] + xy[i1,1] + xy[i2,1]) / 3
            dux_dx_h, dux_dy_h, area = _TriElement._grad_p1(xy, ux, tri)
            duy_dx_h, duy_dy_h, _   = _TriElement._grad_p1(xy, uy, tri)
            err2 += (
                (dux_dx_h - dux_dx(xc, yc, L))**2 +
                (dux_dy_h - dux_dy(xc, yc, L))**2 +
                (duy_dx_h - duy_dx(xc, yc, L))**2 +
                (duy_dy_h - duy_dy(xc, yc, L))**2
            ) * area
        return np.sqrt(err2)

# ==================== SIMULATION ====================
def run_simulation(elem, L, E, nu, nx, ny, dim="2d"):
    root = Sofa.Core.Node("root")
    dofs, nodes_2d, conn = elem.create_sofa_scene(
        root, L=L, E=E, nu=nu, nx=nx, ny=ny, with_visual=False, dim=dim
    )
    Sofa.Simulation.init(root)
    pos0 = dofs.position.array().copy()
    Sofa.Simulation.animate(root, root.dt.value)
    pos1 = dofs.position.array().copy()
    Sofa.Simulation.unload(root)
    ux = pos1[:, 0] - pos0[:, 0]
    uy = pos1[:, 1] - pos0[:, 1]
    return nodes_2d, conn, ux, uy

# ==================== POST-TRAITEMENT ====================
def simulation_ponctuelle(elem, L, E, nu, nx, ny, dim="2d", results_dir=RESULTS_DIR):
    os.makedirs(results_dir, exist_ok=True)
    nodes_2d, conn, ux, uy = run_simulation(elem, L, E, nu, nx, ny, dim=dim)
    l2 = elem.compute_l2(nodes_2d, ux, uy, L, nx, ny, conn)
    h1 = elem.compute_h1(nodes_2d, ux, uy, L, conn)

    xy = nodes_2d[:, :2]
    ux_ref = ux_mms(xy[:, 0], xy[:, 1], L)
    uy_ref = uy_mms(xy[:, 0], xy[:, 1], L)

    hyp = "plane strain" if dim == "3d" else "plane stress"
    print(f"\n[{elem.LABEL}]  {MMS_TYPE.upper()}  dim={dim} ({hyp})  nu={nu}  nx={nx}  L={L}")
    print(f"  L2             = {l2:.4e}")
    print(f"  H1 semi-norm   = {h1:.4e}")
    print(f"  max|ux-ux_mms| = {np.max(np.abs(ux - ux_ref)):.4e}")
    print(f"  max|uy-uy_mms| = {np.max(np.abs(uy - uy_ref)):.4e}")

    mid_j = (ny - 1) // 2
    sl = slice(mid_j * nx, mid_j * nx + nx)
    yc = xy[mid_j * nx, 1]
    xf = np.linspace(0, L, 300)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    for ax, u_sofa, u_fn, lbl, fmt in zip(
        axes,
        [ux[sl], uy[sl]],
        [lambda x: ux_mms(x, yc, L), lambda x: uy_mms(x, yc, L)],
        [r"$u_x$", r"$u_y$"], ["o-", "s-"],
    ):
        ax.plot(xy[sl, 0], u_sofa, fmt, color="tab:orange",
                label=f"SOFA {elem.LABEL} [{dim}]", ms=5)
        ax.plot(xf, u_fn(xf), "--", color="tab:green", label="MMS exact")
        ax.set_xlabel("x")
        ax.set_ylabel(lbl)
        ax.legend()
        ax.grid(True, alpha=0.3)
    plt.suptitle(
        f"MMS {MMS_TYPE.upper()} — {elem.LABEL} [{dim} / {hyp}]  nu={nu}  nx={nx}"
        f"  |L2={l2:.2e}  H1={h1:.2e}"
    )
    plt.tight_layout()
    tag = elem.LABEL.replace(" ", "_")
    plt.savefig(os.path.join(results_dir,
                             f"solution_{tag}_{dim}_nu{nu}_nx{nx}.png"), dpi=150)
    plt.close()

    tris_plot = elem.to_triangles(conn)
    x, y = xy[:, 0], xy[:, 1]
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    for ax, data, title, cmap in [
        (axes[0,0], ux,                r"$u_x$ SOFA",           "gray"),
        (axes[0,1], ux_ref,            r"$u_x$ MMS",            "gray"),
        (axes[0,2], np.abs(ux-ux_ref), r"$|u_x - u_x^{MMS}|$", "gray_r"),
        (axes[1,0], uy,                r"$u_y$ SOFA",           "gray"),
        (axes[1,1], uy_ref,            r"$u_y$ MMS",            "gray"),
        (axes[1,2], np.abs(uy-uy_ref), r"$|u_y - u_y^{MMS}|$", "gray_r"),
    ]:
        tc = ax.tricontourf(x, y, tris_plot.tolist(), data, levels=20, cmap=cmap)
        ax.triplot(x, y, tris_plot.tolist(), "k-", lw=0.3, alpha=0.4)
        plt.colorbar(tc, ax=ax, shrink=0.8)
        ax.set_title(title)
        ax.set_aspect("equal")
        ax.set_xlabel("x"); ax.set_ylabel("y")
    plt.suptitle(
        f"Champs 2D — {elem.LABEL} [{dim} / {hyp}]  {MMS_TYPE.upper()}  nu={nu}  nx={nx}"
    )
    plt.tight_layout()
    plt.savefig(os.path.join(results_dir,
                             f"champs2D_{tag}_{dim}_nu{nu}_nx{nx}.png"), dpi=150)
    plt.close()

# ====================  CONVERGENCE study ====================
def convergence_study(elem_specs, L, E, nu, nx_values, dim="2d", results_dir=RESULTS_DIR):
    os.makedirs(results_dir, exist_ok=True)
    hyp = "plane strain" if dim == "3d" else "plane stress"
    fig_l2, ax_l2 = plt.subplots(figsize=(7, 5))
    fig_h1, ax_h1 = plt.subplots(figsize=(7, 5))

    hs_last = l2_last = h1_last = None

    for spec in elem_specs:
        elem   = spec["elem"]
        label  = spec["label"]
        marker = spec.get("marker", "o")
        color  = spec.get("color", None)

        hs, errs_l2, errs_h1 = [], [], []
        hdr = (f"{'nx':>5} | {'h':>8} | {'L2':>14} | {'ord_L2':>7}"
               f" | {'H1':>14} | {'ord_H1':>7}")
        tag      = label.replace(" ", "_")
        txt_path = os.path.join(results_dir,
                                f"convergence_{tag}_{dim}_nu{nu}.txt")

        print(f"\n── Convergence {label} [{dim} / {hyp}] {MMS_TYPE.upper()} nu={nu} ──\n{hdr}")
        with open(txt_path, "w", encoding="utf-8") as f:
            f.write(f"Convergence MMS {MMS_TYPE.upper()} — {label} [{dim} / {hyp}]  nu={nu}\n{hdr}\n")
            for k, nx in enumerate(nx_values):
                ny = nx
                h  = L / (nx - 1)
                nodes_2d, conn, ux, uy = run_simulation(
                    elem, L, E, nu, nx, ny, dim=dim
                )
                l2 = elem.compute_l2(nodes_2d, ux, uy, L, nx, ny, conn)
                h1 = elem.compute_h1(nodes_2d, ux, uy, L, conn)
                hs.append(h); errs_l2.append(l2); errs_h1.append(h1)

                ord_l2 = (f"{np.log(l2/errs_l2[k-1])/np.log(h/hs[k-1]):.2f}"
                          if k > 0 else "   —  ")
                ord_h1 = (f"{np.log(h1/errs_h1[k-1])/np.log(h/hs[k-1]):.2f}"
                          if k > 0 else "   —  ")
                line = (f"{nx:5d} | {h:8.4f} | {l2:14.6e} | {ord_l2:>7}"
                        f" | {h1:14.6e} | {ord_h1:>7}")
                print(line)
                f.write(line + "\n")

        hs_a = np.array(hs)
        l2_a = np.array(errs_l2)
        h1_a = np.array(errs_h1)
        kw = dict(lw=2, ms=7, color=color)
        ax_l2.loglog(hs_a, l2_a, f"{marker}-",
                     label=f"{label} [{dim}] {MMS_TYPE.upper()} nu={nu}", **kw)
        ax_h1.loglog(hs_a, h1_a, f"{marker}--",
                     label=f"{label} [{dim}] {MMS_TYPE.upper()} nu={nu}", **kw)
        hs_last, l2_last, h1_last = hs_a, l2_a, h1_a

    if hs_last is not None:
        h_ref = np.array([hs_last[0], hs_last[-1]])
        ax_l2.loglog(h_ref, l2_last[0]*(h_ref/hs_last[0])**2,
                     ":", color="gray", lw=1.2, label="O(h²)")
        ax_h1.loglog(h_ref, h1_last[0]*(h_ref/hs_last[0])**1,
                     ":", color="gray", lw=1.2, label="O(h¹)")

    for ax, ylabel, title in [
        (ax_l2, "L2 error",     f"Convergence L2 — MMS {MMS_TYPE.upper()} [{dim} / {hyp}]"),
        (ax_h1, "H1 semi-norm", f"Convergence H1 — MMS {MMS_TYPE.upper()} [{dim} / {hyp}]"),
    ]:
        ax.set_xlabel("h")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend()
        ax.grid(True, alpha=0.3, which="both")

    fig_l2.tight_layout(); fig_h1.tight_layout()
    fig_l2.savefig(os.path.join(results_dir, f"convergence_L2_{dim}_nu{nu}.png"), dpi=150)
    fig_h1.savefig(os.path.join(results_dir, f"convergence_H1_{dim}_nu{nu}.png"), dpi=150)
    plt.close(fig_l2); plt.close(fig_h1)

# ============ MAIN =========================
if __name__ == "__main__":

    print(f"  MMS TYPE : {MMS_TYPE.upper()}")

    L  = 1.0
    E  = 1e6
    nx_values = [10, 20, 50, 80, 100, 120, 150, 200, 250 ]

    element_quad = _QuadElement()
    element_tri  = _TriElement()

    specs = [
        {"elem": element_quad, "label": "Q1 quad", "marker": "o", "color": "C0"},
        {"elem": element_tri,  "label": "P1 tri",  "marker": "s", "color": "C1"},
    ]

    # =========== Plane stress vec2d =======================
    DIM = "2d"
    for nu in [0.0, 0.3, 0.49]:
        print(f"\n  PoissonRatio = {nu} ")
        simulation_ponctuelle(element_quad, L, E, nu, nx=20, ny=20, dim=DIM)
        simulation_ponctuelle(element_tri,  L, E, nu, nx=20, ny=20, dim=DIM)
        convergence_study(specs, L, E, nu, nx_values, dim=DIM)

    # ============== Plane strain vec3d  ======================
    DIM = "3d"
    
    for nu in [0.0, 0.3, 0.49]:
        print(f"\n  PoissonRatio = {nu}  ")
        simulation_ponctuelle(element_quad, L, E, nu, nx=20, ny=20, dim=DIM)
        simulation_ponctuelle(element_tri,  L, E, nu, nx=20, ny=20, dim=DIM)
        convergence_study(specs, L, E, nu, nx_values, dim=DIM)
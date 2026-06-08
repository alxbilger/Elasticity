import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import logging

import Sofa
import Sofa.Core
import Sofa.Simulation
import SofaRuntime
import gmsh

from gmsh_generate_beam2D import generate_beam2D

"""
Incompressible ===> div u(x,y) = 0,
MMS realistic : u_x = sin(pi*x)cos(pi*y);   u_y = -cos(pi*x)sin(pi*y)

Locking ratio (Bathe 1996 / Zienkiewicz & Taylor) :
    LR(nu, h) = |u_h|_{H1} / |u_exact|_{H1}
Un élément verrouillé -> LR -> 0 quand nu -> 0.5.
Un élément libre du locking -> LR -> 1 quelle que soit nu.
"""

RESULTS_DIR = "results_incompressible_q1p1"
os.makedirs(RESULTS_DIR, exist_ok=True)

# ===== Silence SOFA console output =====
logging.getLogger("Sofa").setLevel(logging.ERROR)

# ===== Gauss points & weights =====
GAUSS_PTS = np.array([-1.0 / np.sqrt(3), 1.0 / np.sqrt(3)])
GAUSS_WTS = np.array([1.0, 1.0])

# ===== Triangle Gauss rule (3-point, ordre 2) =====
_TRI_BARY = np.array([
    [2/3, 1/6, 1/6],
    [1/6, 2/3, 1/6],
    [1/6, 1/6, 2/3],
])


# ======= Lamé coefficients =========================
def lame(E, nu, dim="2d"):
    if dim == "3d" and np.isclose(nu, 0.5):
        raise ValueError(
            "nu=0.5 in plane strain ===> division by zero. Use nu < 0.5"
        )
    if dim == "2d":
        lam = E * nu / (1 - nu**2)
    else:
        lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    return lam, mu


# ==================== MMS & derivatives ==========================
def ux_mms_incomp(x, y, L):  return  np.sin(np.pi*x/L) * np.cos(np.pi*y/L)
def uy_mms_incomp(x, y, L):  return -np.cos(np.pi*x/L) * np.sin(np.pi*y/L)

def dux_dx_incomp(x, y, L):  return  (np.pi/L)*np.cos(np.pi*x/L)*np.cos(np.pi*y/L)
def dux_dy_incomp(x, y, L):  return -(np.pi/L)*np.sin(np.pi*x/L)*np.sin(np.pi*y/L)
def duy_dx_incomp(x, y, L):  return  (np.pi/L)*np.sin(np.pi*x/L)*np.sin(np.pi*y/L)
def duy_dy_incomp(x, y, L):  return -(np.pi/L)*np.cos(np.pi*x/L)*np.cos(np.pi*y/L)

def fx_body_incomp(x, y, E, nu, L, dim="2d"):
    return (np.pi**2 * E / (L**2 * (nu + 1))) * np.sin(np.pi*x/L) * np.cos(np.pi*y/L)

def fy_body_incomp(x, y, E, nu, L, dim="2d"):
    return -(np.pi**2 * E / (L**2 * (nu + 1))) * np.sin(np.pi*y/L) * np.cos(np.pi*x/L)

def sigma_xx_incomp(x, y, E, nu, L, dim="2d"):
    _, mu = lame(E, nu, dim)
    return 2*mu * dux_dx_incomp(x, y, L)

def sigma_yy_incomp(x, y, E, nu, L, dim="2d"):
    _, mu = lame(E, nu, dim)
    return 2*mu * duy_dy_incomp(x, y, L)

def sigma_xy_incomp(x, y, E, nu, L, dim="2d"):
    _, mu = lame(E, nu, dim)
    return mu * (dux_dy_incomp(x, y, L) + duy_dx_incomp(x, y, L))

def traction_incomp(x, y, nx_c, ny_c, E, nu, L, dim="2d"):
    sxx = sigma_xx_incomp(x, y, E, nu, L, dim)
    syy = sigma_yy_incomp(x, y, E, nu, L, dim)
    sxy = sigma_xy_incomp(x, y, E, nu, L, dim)
    return sxx*nx_c + sxy*ny_c, sxy*nx_c + syy*ny_c


# =================== Norme H1 exacte de la solution MMS ====================

def compute_h1_exact_tri(nodes_2d, conn, L):
    """
    ||u_exact||_{H1}  calculée par quadrature sur les triangles du maillage.
    Utilisée au dénominateur du locking ratio.
    """
    xy   = nodes_2d[:, :2]
    val2 = 0.0
    for tri in conn:
        i0, i1, i2 = tri
        x0,y0 = xy[i0]; x1,y1 = xy[i1]; x2,y2 = xy[i2]
        area  = abs((x1-x0)*(y2-y0) - (x2-x0)*(y1-y0)) / 2.0
        xc    = (x0 + x1 + x2) / 3.0
        yc    = (y0 + y1 + y2) / 3.0
        val2 += (
            dux_dx_incomp(xc, yc, L)**2 +
            dux_dy_incomp(xc, yc, L)**2 +
            duy_dx_incomp(xc, yc, L)**2 +
            duy_dy_incomp(xc, yc, L)**2
        ) * area
    return np.sqrt(val2)

def compute_h1_exact_quad(nodes_2d, conn, L):
    """
    ||u_exact||_{H1}  calculée sur les quads (2x2 Gauss).
    """
    xy   = nodes_2d[:, :2]
    val2 = 0.0
    for quad in conn:
        xe, ye = xy[quad, 0], xy[quad, 1]
        for xi, wi in zip(GAUSS_PTS, GAUSS_WTS):
            for eta, wj in zip(GAUSS_PTS, GAUSS_WTS):
                N, dN_dxi, dN_deta = _QuadElement._shape(xi, eta)
                J    = np.array([[dN_dxi@xe, dN_dxi@ye],
                                 [dN_deta@xe, dN_deta@ye]])
                detJ = np.linalg.det(J)
                Jinv = np.linalg.inv(J)
                xg, yg = N@xe, N@ye
                val2 += (
                    dux_dx_incomp(xg, yg, L)**2 +
                    dux_dy_incomp(xg, yg, L)**2 +
                    duy_dx_incomp(xg, yg, L)**2 +
                    duy_dy_incomp(xg, yg, L)**2
                ) * wi * wj * detJ
    return np.sqrt(val2)


# =================== Mesh helpers ============================================

def _sofa_template(dim):
    return "Vec3d" if dim == "3d" else "Vec2d"

def _add_dirichlet(Beam, nodes_2d, L, dim, with_visual):
    tmpl = _sofa_template(dim)
    xy   = nodes_2d[:, :2]
    eps  = 1e-10

    idx_corner = [k for k, (xk, yk) in enumerate(xy)
                  if xk < eps and yk < eps]
    idx_rot    = [k for k, (xk, yk) in enumerate(xy)
                  if xk > L - eps and yk < eps]

    if dim == "2d":
        Beam.addObject("PartialFixedProjectiveConstraint",
                       name="fix_corner", template=tmpl,
                       indices=" ".join(map(str, idx_corner)),
                       fixedDirections="1 1")
        Beam.addObject("PartialFixedProjectiveConstraint",
                       name="fix_rot", template=tmpl,
                       indices=" ".join(map(str, idx_rot)),
                       fixedDirections="0 1")
    else:
        all_idx = " ".join(str(k) for k in range(len(nodes_2d)))
        Beam.addObject("PartialFixedProjectiveConstraint",
                       name="fix_corner", template=tmpl,
                       indices=" ".join(map(str, idx_corner)),
                       fixedDirections="1 1 0")
        Beam.addObject("PartialFixedProjectiveConstraint",
                       name="fix_rot", template=tmpl,
                       indices=" ".join(map(str, idx_rot)),
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


# =================== Gmsh mesh reader ========================================

def load_gmsh_mesh(msh_path, dim="2d"):
    """
    Lit un fichier .msh (MSH v2) et retourne :
      - nodes_2d      : (N, 2|3)
      - tris          : (T, 3)  indices 0-based
      - boundary_edges: dict {left/right/bottom/top -> [(k0,k1),...]}
    Physical groups : left=1, bottom=2, right=3, top=4
    """
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)   # silence gmsh stdout
    gmsh.open(msh_path)

    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    xyz        = node_coords.reshape(-1, 3)
    tag_to_idx = {int(t): i for i, t in enumerate(node_tags)}
    n_nodes    = len(node_tags)

    if dim == "3d":
        nodes_2d = np.column_stack([xyz[:, 0], xyz[:, 1], np.zeros(n_nodes)])
    else:
        nodes_2d = xyz[:, :2].copy()

    elem_types, _, elem_node_tags = gmsh.model.mesh.getElements(dim=2)
    tris = []
    for etype, entags in zip(elem_types, elem_node_tags):
        if etype == 2:  # tri3
            conn = entags.reshape(-1, 3)
            for tri in conn:
                tris.append([tag_to_idx[int(tri[0])],
                              tag_to_idx[int(tri[1])],
                              tag_to_idx[int(tri[2])]])
    tris = np.array(tris, dtype=int)

    label_map = {1: "left", 2: "bottom", 3: "right", 4: "top"}
    boundary_edges = {name: [] for name in label_map.values()}
    for phys_tag, name in label_map.items():
        try:
            entities = gmsh.model.getEntitiesForPhysicalGroup(1, phys_tag)
        except Exception:
            continue
        for ent in entities:
            etypes, _, entags_1d = gmsh.model.mesh.getElements(dim=1, tag=int(ent))
            for etype, enodes in zip(etypes, entags_1d):
                if etype == 1:
                    edges = enodes.reshape(-1, 2)
                    for e in edges:
                        boundary_edges[name].append(
                            (tag_to_idx[int(e[0])], tag_to_idx[int(e[1])])
                        )
    gmsh.finalize()
    return nodes_2d, tris, boundary_edges


# ========================   Q1 Quad element ==================================

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
    def get_nodes_2d(L, nx, ny, dim="2d"):
        dx, dy = L / (nx - 1), L / (ny - 1)
        if dim == "3d":
            pts = [[i*dx, j*dy, 0.0] for j in range(ny) for i in range(nx)]
        else:
            pts = [[i*dx, j*dy]      for j in range(ny) for i in range(nx)]
        return np.array(pts, dtype=float)

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
                    J    = np.array([[dN_dxi@xe, dN_dxi@ye],
                                     [dN_deta@xe, dN_deta@ye]])
                    detJ = np.linalg.det(J)
                    xg, yg = N@xe, N@ye
                    fx = fx_body_incomp(xg, yg, E, nu, L, dim)
                    fy = fy_body_incomp(xg, yg, E, nu, L, dim)
                    w  = wi * wj * detJ
                    for a, node in enumerate(quad):
                        F[node, 0] += N[a] * fx * w
                        F[node, 1] += N[a] * fy * w

        def _edge_gauss(k0, k1, nx_c, ny_c):
            x0, y0 = xy[k0]; x1, y1 = xy[k1]
            Le = np.hypot(x1-x0, y1-y0)
            for xi, wi in zip(GAUSS_PTS, GAUSS_WTS):
                t      = 0.5*(xi + 1)
                xg, yg = (1-t)*x0 + t*x1, (1-t)*y0 + t*y1
                Tx, Ty = traction_incomp(xg, yg, nx_c, ny_c, E, nu, L, dim)
                N0, N1 = 0.5*(1-xi), 0.5*(1+xi)
                w      = wi * Le / 2.0
                F[k0, 0] += N0*Tx*w;  F[k0, 1] += N0*Ty*w
                F[k1, 0] += N1*Tx*w;  F[k1, 1] += N1*Ty*w

        for j in range(ny - 1):
            _edge_gauss(j*nx, (j+1)*nx, -1., 0.)
        for j in range(ny - 1):
            _edge_gauss(j*nx+(nx-1), (j+1)*nx+(nx-1), +1., 0.)
        for i in range(nx - 1):
            _edge_gauss(i, i+1, 0., -1.)
        for i in range(nx - 1):
            _edge_gauss((ny-1)*nx+i, (ny-1)*nx+i+1, 0., +1.)

        return F

    @staticmethod
    def create_sofa_scene(rootNode, L=1.0, E=1e6, nu=0.3,
                          nx=10, ny=10, with_visual=True, dim="2d"):
        tmpl = _sofa_template(dim)
        rootNode.addObject("RequiredPlugin", pluginName=[
            "Sofa.Component.Visual",
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

        nodes_2d = _QuadElement.get_nodes_2d(L, nx, ny, dim=dim)
        quads_np = _QuadElement.get_connectivity(nx, ny)

        Beam = rootNode.addChild("Beam")
        Beam.addObject("NewtonRaphsonSolver", name="newtonSolver",
                       maxNbIterationsNewton=10,
                       absoluteResidualStoppingThreshold=1e-10,
                       printLog=False)
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
                    J    = np.array([[dN_dxi@xe,  dN_dxi@ye],
                                     [dN_deta@xe, dN_deta@ye]])
                    detJ = np.linalg.det(J)
                    xg, yg = N@xe, N@ye
                    err2 += (
                        (N@ux[quad] - ux_mms_incomp(xg, yg, L))**2
                      + (N@uy[quad] - uy_mms_incomp(xg, yg, L))**2
                    ) * wi * wj * detJ
        return np.sqrt(err2)

    @staticmethod
    def compute_h1(nodes_2d, ux, uy, L, conn, **_):
        xy, err2 = nodes_2d[:, :2], 0.0
        for quad in conn:
            xe, ye = xy[quad, 0], xy[quad, 1]
            for xi, wi in zip(GAUSS_PTS, GAUSS_WTS):
                for eta, wj in zip(GAUSS_PTS, GAUSS_WTS):
                    N, dN_dxi, dN_deta = _QuadElement._shape(xi, eta)
                    J    = np.array([[dN_dxi@xe,  dN_dxi@ye],
                                     [dN_deta@xe, dN_deta@ye]])
                    detJ = np.linalg.det(J)
                    Jinv = np.linalg.inv(J)
                    xg, yg = N@xe, N@ye
                    dN_dx = Jinv[0,0]*dN_dxi + Jinv[1,0]*dN_deta
                    dN_dy = Jinv[0,1]*dN_dxi + Jinv[1,1]*dN_deta
                    err2 += (
                        (dN_dx@ux[quad] - dux_dx_incomp(xg, yg, L))**2 +
                        (dN_dy@ux[quad] - dux_dy_incomp(xg, yg, L))**2 +
                        (dN_dx@uy[quad] - duy_dx_incomp(xg, yg, L))**2 +
                        (dN_dy@uy[quad] - duy_dy_incomp(xg, yg, L))**2
                    ) * wi * wj * detJ
        return np.sqrt(err2)

    @staticmethod
    def compute_h1_sofa(nodes_2d, ux, uy, conn):
        """
        ||u_h||_{H1}  (semi-norme) — numérateur du locking ratio pour Q1.
        """
        xy, val2 = nodes_2d[:, :2], 0.0
        for quad in conn:
            xe, ye = xy[quad, 0], xy[quad, 1]
            for xi, wi in zip(GAUSS_PTS, GAUSS_WTS):
                for eta, wj in zip(GAUSS_PTS, GAUSS_WTS):
                    N, dN_dxi, dN_deta = _QuadElement._shape(xi, eta)
                    J    = np.array([[dN_dxi@xe,  dN_dxi@ye],
                                     [dN_deta@xe, dN_deta@ye]])
                    detJ = np.linalg.det(J)
                    Jinv = np.linalg.inv(J)
                    dN_dx = Jinv[0,0]*dN_dxi + Jinv[1,0]*dN_deta
                    dN_dy = Jinv[0,1]*dN_dxi + Jinv[1,1]*dN_deta
                    val2 += (
                        (dN_dx@ux[quad])**2 + (dN_dy@ux[quad])**2 +
                        (dN_dx@uy[quad])**2 + (dN_dy@uy[quad])**2
                    ) * wi * wj * detJ
        return np.sqrt(val2)


# ====================  P1 Triangle element — maillage via Gmsh ===============

class _TriElement:
    LABEL = "P1 tri"

    @staticmethod
    def load_mesh(L, nx, ny, dim="2d"):
        filename = f"mesh_tri_nx{nx}_ny{ny}.msh"
        msh_path, _, _ = generate_beam2D(
            length=L, height=L, nx=nx, ny=ny, filename=filename,
        )
        nodes_2d, tris, boundary_edges = load_gmsh_mesh(msh_path, dim=dim)
        return nodes_2d, tris, boundary_edges

    @staticmethod
    def _grad_p1(xy, u, tri):
        i0, i1, i2 = tri
        x0,y0 = xy[i0];  x1,y1 = xy[i1];  x2,y2 = xy[i2]
        A2   = (x1-x0)*(y2-y0) - (x2-x0)*(y1-y0)
        dudx = (u[i0]*(y1-y2) + u[i1]*(y2-y0) + u[i2]*(y0-y1)) / A2
        dudy = (u[i0]*(x2-x1) + u[i1]*(x0-x2) + u[i2]*(x1-x0)) / A2
        return dudx, dudy, abs(A2) / 2.0

    @staticmethod
    def to_triangles(conn):
        return conn

    @staticmethod
    def compute_nodal_forces(nodes_2d, L, E, nu, tris, boundary_edges, dim="2d"):
        xy = nodes_2d[:, :2]
        F  = np.zeros((len(xy), 2))

        for tri in tris:
            i0, i1, i2 = tri
            x0,y0 = xy[i0]; x1,y1 = xy[i1]; x2,y2 = xy[i2]
            xc   = (x0 + x1 + x2) / 3.0
            yc   = (y0 + y1 + y2) / 3.0
            area = abs((x1-x0)*(y2-y0) - (x2-x0)*(y1-y0)) / 2.0
            fx   = fx_body_incomp(xc, yc, E, nu, L, dim)
            fy   = fy_body_incomp(xc, yc, E, nu, L, dim)
            for node in [i0, i1, i2]:
                F[node, 0] += fx * area / 3.0
                F[node, 1] += fy * area / 3.0

        normal_map = {
            "left":   (-1.,  0.), "right":  (+1.,  0.),
            "bottom": ( 0., -1.), "top":    ( 0., +1.),
        }

        def _edge_gauss(k0, k1, nx_c, ny_c):
            x0, y0 = xy[k0]; x1, y1 = xy[k1]
            Le = np.hypot(x1-x0, y1-y0)
            for xi, wi in zip(GAUSS_PTS, GAUSS_WTS):
                t      = 0.5*(xi + 1)
                xg, yg = (1-t)*x0 + t*x1, (1-t)*y0 + t*y1
                Tx, Ty = traction_incomp(xg, yg, nx_c, ny_c, E, nu, L, dim)
                N0, N1 = 0.5*(1-xi), 0.5*(1+xi)
                w      = wi * Le / 2.0
                F[k0, 0] += N0*Tx*w;  F[k0, 1] += N0*Ty*w
                F[k1, 0] += N1*Tx*w;  F[k1, 1] += N1*Ty*w

        for bname, (nx_c, ny_c) in normal_map.items():
            for (k0, k1) in boundary_edges[bname]:
                _edge_gauss(k0, k1, nx_c, ny_c)

        return F

    @staticmethod
    def create_sofa_scene(rootNode, L=1.0, E=1e6, nu=0.3,
                          nx=10, ny=10, with_visual=True, dim="2d"):
        tmpl = _sofa_template(dim)
        rootNode.addObject("RequiredPlugin", pluginName=[
            "Sofa.Component.Visual",
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

        nodes_2d, tris_np, boundary_edges = _TriElement.load_mesh(L, nx, ny, dim=dim)

        Beam = rootNode.addChild("Beam")
        Beam.addObject("NewtonRaphsonSolver", name="newtonSolver",
                       maxNbIterationsNewton=10,
                       absoluteResidualStoppingThreshold=1e-10,
                       printLog=False)
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
        F_xy = _TriElement.compute_nodal_forces(
            nodes_2d, L, E, nu, tris_np, boundary_edges, dim
        )
        idx, frc = _pack_forces(F_xy, dim, len(nodes_2d))
        Beam.addObject("ConstantForceField", name="MMS_forces", template=tmpl,
                       indices=idx, forces=frc)
        return dofs, nodes_2d, tris_np

    @staticmethod
    def compute_l2(nodes_2d, ux, uy, L, nx, ny, conn):
        xy   = nodes_2d[:, :2]
        err2 = 0.0
        for tri in conn:
            i0, i1, i2 = tri
            xy_e = xy[[i0, i1, i2]]
            ux_e = ux[[i0, i1, i2]]
            uy_e = uy[[i0, i1, i2]]
            _, _, area = _TriElement._grad_p1(xy, ux, tri)
            for lam in _TRI_BARY:
                xg  = lam @ xy_e[:, 0]
                yg  = lam @ xy_e[:, 1]
                uxg = lam @ ux_e
                uyg = lam @ uy_e
                err2 += (
                    (uxg - ux_mms_incomp(xg, yg, L))**2
                  + (uyg - uy_mms_incomp(xg, yg, L))**2
                ) * area / 3.0
        return np.sqrt(err2)

    @staticmethod
    def compute_h1(nodes_2d, ux, uy, L, conn, **_):
        xy   = nodes_2d[:, :2]
        err2 = 0.0
        for tri in conn:
            i0, i1, i2 = tri
            xc = (xy[i0,0] + xy[i1,0] + xy[i2,0]) / 3.0
            yc = (xy[i0,1] + xy[i1,1] + xy[i2,1]) / 3.0
            dux_dx_h, dux_dy_h, area = _TriElement._grad_p1(xy, ux, tri)
            duy_dx_h, duy_dy_h, _   = _TriElement._grad_p1(xy, uy, tri)
            err2 += (
                (dux_dx_h - dux_dx_incomp(xc, yc, L))**2 +
                (dux_dy_h - dux_dy_incomp(xc, yc, L))**2 +
                (duy_dx_h - duy_dx_incomp(xc, yc, L))**2 +
                (duy_dy_h - duy_dy_incomp(xc, yc, L))**2
            ) * area
        return np.sqrt(err2)

    @staticmethod
    def compute_h1_sofa(nodes_2d, ux, uy, conn):
        """
        ||u_h||_{H1}  (semi-norme) — numérateur du locking ratio pour P1.
        """
        xy   = nodes_2d[:, :2]
        val2 = 0.0
        for tri in conn:
            dux_dx_h, dux_dy_h, area = _TriElement._grad_p1(xy, ux, tri)
            duy_dx_h, duy_dy_h, _   = _TriElement._grad_p1(xy, uy, tri)
            val2 += (
                dux_dx_h**2 + dux_dy_h**2 +
                duy_dx_h**2 + duy_dy_h**2
            ) * area
        return np.sqrt(val2)


# ============================== SIMU =========================================

def run_simulation(elem, L, E, nu, nx, ny, dim="2d"):
    """Lance une simulation SOFA et retourne (nodes_2d, conn, ux, uy)."""
    # Redirect stdout/stderr to suppress SOFA console spam
    devnull = open(os.devnull, 'w')
    old_stdout, old_stderr = sys.stdout, sys.stderr
    sys.stdout = sys.stderr = devnull
    try:
        root = Sofa.Core.Node("root")
        dofs, nodes_2d, conn = elem.create_sofa_scene(
            root, L=L, E=E, nu=nu, nx=nx, ny=ny, with_visual=False, dim=dim
        )
        Sofa.Simulation.init(root)
        pos0 = dofs.position.array().copy()
        Sofa.Simulation.animate(root, root.dt.value)
        pos1 = dofs.position.array().copy()
        Sofa.Simulation.unload(root)
    finally:
        sys.stdout, sys.stderr = old_stdout, old_stderr
        devnull.close()
    ux = pos1[:, 0] - pos0[:, 0]
    uy = pos1[:, 1] - pos0[:, 1]
    return nodes_2d, conn, ux, uy


def simulation_ponctuelle(elem, L, E, nu, nx, ny, dim="2d", results_dir=RESULTS_DIR):
    os.makedirs(results_dir, exist_ok=True)
    nodes_2d, conn, ux, uy = run_simulation(elem, L, E, nu, nx, ny, dim=dim)
    l2 = elem.compute_l2(nodes_2d, ux, uy, L, nx, ny, conn)
    h1 = elem.compute_h1(nodes_2d, ux, uy, L, conn)

    # --- locking ratio ---
    h1_sofa = elem.compute_h1_sofa(nodes_2d, ux, uy, conn)
    if isinstance(elem, _TriElement):
        h1_exact = compute_h1_exact_tri(nodes_2d, conn, L)
    else:
        h1_exact = compute_h1_exact_quad(nodes_2d, conn, L)
    locking_ratio = h1_sofa / h1_exact if h1_exact > 0 else float('nan')

    xy     = nodes_2d[:, :2]
    ux_ref = ux_mms_incomp(xy[:, 0], xy[:, 1], L)
    uy_ref = uy_mms_incomp(xy[:, 0], xy[:, 1], L)

    hyp = "plane strain" if dim == "3d" else "plane stress"
    print(f"  L2             = {l2:.4e}")
    print(f"  H1 error       = {h1:.4e}")
    print(f"  Locking ratio  = {locking_ratio:.4f}  "
          f"(1=no locking, 0=full locking)")
    print(f"  max|ux-ux_mms| = {np.max(np.abs(ux - ux_ref)):.4e}")
    print(f"  max|uy-uy_mms| = {np.max(np.abs(uy - uy_ref)):.4e}")

    yc      = L / 2.0
    eps_cut = L / (2 * (ny - 1)) if ny > 1 else 0.05
    mask    = np.abs(xy[:, 1] - yc) < eps_cut
    xf      = np.linspace(0, L, 300)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    for ax, u_sofa, u_fn, lbl in zip(
        axes,
        [ux[mask], uy[mask]],
        [lambda x: ux_mms_incomp(x, yc, L), lambda x: uy_mms_incomp(x, yc, L)],
        [r"$u_x$", r"$u_y$"],
    ):
        ax.scatter(xy[mask, 0], u_sofa, c="tab:orange",
                   label=f"SOFA {elem.LABEL} [{dim}]", s=20)
        ax.plot(xf, u_fn(xf), "--", color="tab:green", label="MMS exact")
        ax.set_xlabel("x"); ax.set_ylabel(lbl)
        ax.legend(); ax.grid(True, alpha=0.3)
    plt.suptitle(
        f"MMS INCOMPRESSIBLE — {elem.LABEL} [{dim}/{hyp}]  nu={nu}  nx={nx}"
        f"  |L2={l2:.2e}  H1={h1:.2e}  LR={locking_ratio:.3f}"
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
        (axes[0,0], ux,                r"$u_x$ SOFA",          "gray"),
        (axes[0,1], ux_ref,            r"$u_x$ MMS",           "gray"),
        (axes[0,2], np.abs(ux-ux_ref), r"$|u_x-u_x^{MMS}|$",  "gray_r"),
        (axes[1,0], uy,                r"$u_y$ SOFA",          "gray"),
        (axes[1,1], uy_ref,            r"$u_y$ MMS",           "gray"),
        (axes[1,2], np.abs(uy-uy_ref), r"$|u_y-u_y^{MMS}|$",  "gray_r"),
    ]:
        tc = ax.tricontourf(x, y, tris_plot.tolist(), data, levels=20, cmap=cmap)
        ax.triplot(x, y, tris_plot.tolist(), "k-", lw=0.3, alpha=0.4)
        plt.colorbar(tc, ax=ax, shrink=0.8)
        ax.set_title(title); ax.set_aspect("equal")
        ax.set_xlabel("x"); ax.set_ylabel("y")
    plt.suptitle(
        f"Champs 2D — {elem.LABEL} [{dim}/{hyp}]  INCOMPRESSIBLE  "
        f"nu={nu}  nx={nx}  LR={locking_ratio:.3f}"
    )
    plt.tight_layout()
    plt.savefig(os.path.join(results_dir,
                f"champs2D_{tag}_{dim}_nu{nu}_nx{nx}.png"), dpi=150)
    plt.close()


def convergence_study(elem_specs, L, E, nu, nx_values, dim="2d",
                      results_dir=RESULTS_DIR):
    os.makedirs(results_dir, exist_ok=True)
    hyp = "plane strain" if dim == "3d" else "plane stress"
    fig_l2, ax_l2   = plt.subplots(figsize=(7, 5))
    fig_h1, ax_h1   = plt.subplots(figsize=(7, 5))
    fig_lr, ax_lr   = plt.subplots(figsize=(7, 5))   # locking ratio
    hs_last = l2_last = h1_last = None

    for spec in elem_specs:
        elem   = spec["elem"]
        label  = spec["label"]
        marker = spec.get("marker", "o")
        color  = spec.get("color", None)

        hs, errs_l2, errs_h1, lock_ratios = [], [], [], []
        hdr = (f"{'nx':>5} | {'h':>8} | {'L2':>14} | {'ord_L2':>7}"
               f" | {'H1':>14} | {'ord_H1':>7} | {'LR':>8}")
        tag      = label.replace(" ", "_")
        txt_path = os.path.join(results_dir,
                                f"convergence_{tag}_{dim}_nu{nu}.txt")

        print(f"\n── Convergence {label} [{dim}/{hyp}] INCOMP nu={nu} ──\n{hdr}")
        with open(txt_path, "w", encoding="utf-8") as f:
            f.write(
                f"Convergence MMS INCOMPRESSIBLE — {label} [{dim}/{hyp}]  nu={nu}\n"
                f"Locking ratio = ||u_h||_H1 / ||u_exact||_H1\n{hdr}\n"
            )
            for k, nx in enumerate(nx_values):
                ny = nx
                h  = L / (nx - 1)
                nodes_2d, conn, ux, uy = run_simulation(
                    elem, L, E, nu, nx, ny, dim=dim
                )
                l2 = elem.compute_l2(nodes_2d, ux, uy, L, nx, ny, conn)
                h1 = elem.compute_h1(nodes_2d, ux, uy, L, conn)

                # locking ratio
                h1_sofa = elem.compute_h1_sofa(nodes_2d, ux, uy, conn)
                if isinstance(elem, _TriElement):
                    h1_ex = compute_h1_exact_tri(nodes_2d, conn, L)
                else:
                    h1_ex = compute_h1_exact_quad(nodes_2d, conn, L)
                lr = h1_sofa / h1_ex if h1_ex > 0 else float('nan')

                hs.append(h); errs_l2.append(l2)
                errs_h1.append(h1); lock_ratios.append(lr)

                ord_l2 = (f"{np.log(l2/errs_l2[k-1])/np.log(h/hs[k-1]):.2f}"
                          if k > 0 else "   —  ")
                ord_h1 = (f"{np.log(h1/errs_h1[k-1])/np.log(h/hs[k-1]):.2f}"
                          if k > 0 else "   —  ")
                line = (f"{nx:5d} | {h:8.4f} | {l2:14.6e} | {ord_l2:>7}"
                        f" | {h1:14.6e} | {ord_h1:>7} | {lr:8.4f}")
                print(line); f.write(line + "\n")

        hs_a  = np.array(hs)
        l2_a  = np.array(errs_l2)
        h1_a  = np.array(errs_h1)
        lr_a  = np.array(lock_ratios)
        kw = dict(lw=2, ms=7, color=color)
        ax_l2.loglog(hs_a, l2_a, f"{marker}-",
                     label=f"{label} [{dim}] nu={nu}", **kw)
        ax_h1.loglog(hs_a, h1_a, f"{marker}--",
                     label=f"{label} [{dim}] nu={nu}", **kw)
        ax_lr.semilogx(hs_a, lr_a, f"{marker}-",
                       label=f"{label} [{dim}] nu={nu}", **kw)
        hs_last, l2_last, h1_last = hs_a, l2_a, h1_a

    if hs_last is not None:
        h_ref = np.array([hs_last[0], hs_last[-1]])
        ax_l2.loglog(h_ref, l2_last[0]*(h_ref/hs_last[0])**2,
                     ":", color="gray", lw=1.2, label="O(h²)")
        ax_h1.loglog(h_ref, h1_last[0]*(h_ref/hs_last[0])**1,
                     ":", color="gray", lw=1.2, label="O(h¹)")

    for ax, ylabel, title in [
        (ax_l2, "L2 error",     f"Convergence L2 [{dim}/{hyp}] nu={nu}"),
        (ax_h1, "H1 semi-norm", f"Convergence H1 [{dim}/{hyp}] nu={nu}"),
    ]:
        ax.set_xlabel("h"); ax.set_ylabel(ylabel); ax.set_title(title)
        ax.legend(); ax.grid(True, alpha=0.3, which="both")

    ax_lr.axhline(1.0, color="gray", ls=":", lw=1.2, label="LR=1 (no locking)")
    ax_lr.set_xlabel("h")
    ax_lr.set_ylabel(r"Locking ratio  $\|u_h\|_{H^1} / \|u_{ex}\|_{H^1}$")
    ax_lr.set_title(f"Locking ratio [{dim}/{hyp}] nu={nu}")
    ax_lr.set_ylim(0, 1.1)
    ax_lr.legend(); ax_lr.grid(True, alpha=0.3, which="both")

    fig_l2.tight_layout(); fig_h1.tight_layout(); fig_lr.tight_layout()
    fig_l2.savefig(os.path.join(results_dir,
                   f"convergence_L2_{dim}_nu{nu}.png"), dpi=150)
    fig_h1.savefig(os.path.join(results_dir,
                   f"convergence_H1_{dim}_nu{nu}.png"), dpi=150)
    fig_lr.savefig(os.path.join(results_dir,
                   f"locking_ratio_{dim}_nu{nu}.png"), dpi=150)
    plt.close(fig_l2); plt.close(fig_h1); plt.close(fig_lr)


# ============================ Main ===========================================
if __name__ == "__main__":
    L         = 1.0
    E         = 1e6
    nx_values = [10, 20, 50, 80, 100, 120, 150, 200]

    element_quad = _QuadElement()
    element_tri  = _TriElement()

    specs = [
        {"elem": element_quad, "label": "Q1 quad", "marker": "o", "color": "C0"},
        {"elem": element_tri,  "label": "P1 tri",  "marker": "s", "color": "C1"},
    ]

    DIM = "2d"
    for nu in [0.0, 0.3, 0.49, 0.4999]:
        print(f"\n  nu = {nu}  [{DIM}]")
        simulation_ponctuelle(element_quad, L, E, nu, nx=20, ny=20, dim=DIM)
        simulation_ponctuelle(element_tri,  L, E, nu, nx=20, ny=20, dim=DIM)
        convergence_study(specs, L, E, nu, nx_values, dim=DIM)

    DIM = "3d"
    for nu in [0.0, 0.3, 0.49, 0.4999]:
        print(f"\n  nu = {nu}  [{DIM}]")
        simulation_ponctuelle(element_quad, L, E, nu, nx=20, ny=20, dim=DIM)
        simulation_ponctuelle(element_tri,  L, E, nu, nx=20, ny=20, dim=DIM)
        convergence_study(specs, L, E, nu, nx_values, dim=DIM)

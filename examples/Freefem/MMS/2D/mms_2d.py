import numpy as np
import matplotlib.pyplot as plt
import os

import Sofa
import Sofa.Core
import Sofa.Simulation
import SofaRuntime

RESULTS_DIR = "results_mms2d_Q1P1"
os.makedirs(RESULTS_DIR, exist_ok=True)

GAUSS_PTS = np.array([-1.0 / np.sqrt(3), 1.0 / np.sqrt(3)])
GAUSS_WTS = np.array([1.0, 1.0])



#============  MMS manufactured solution & derivatives =====================


def ux_mms(x, y, L): return x**2 * (L - x) / L**2

def uy_mms(x, y, L): return x * (L - x) * y / L**2

def dux_dx(x, y, L): return (2*x*L - 3*x**2) / L**2

def dux_dy(x, y, L): return np.zeros_like(np.asarray(x, float))

def duy_dx(x, y, L): return (L - 2*x) * y / L**2

def duy_dy(x, y, L): return x * (L - x) / L**2

def sigma_xx(x, y, E, L): return E * dux_dx(x, y, L)

def sigma_yy(x, y, E, L): return E * duy_dy(x, y, L)

def sigma_xy(x, y, E, L): return E * 0.5 * duy_dx(x, y, L)

def fx_body(x, y, E, L): return -(E*(2*L - 6*x)/L**2 + E*(L - 2*x)/(2*L**2))

def fy_body(x, y, E, L): return E * y / L**2

def traction(x, y, nx_c, ny_c, E, L):
    sxx = sigma_xx(x, y, E, L)
    syy = sigma_yy(x, y, E, L)
    sxy = sigma_xy(x, y, E, L)
    return sxx*nx_c + sxy*ny_c, sxy*nx_c + syy*ny_c



# ============== Mesh helpers ============================


def get_nodes_2d(L, nx, ny, dim="2d"):
    dx = L / (nx - 1)
    dy = L / (ny - 1)
    if dim == "3d":
        pts = [[i*dx, j*dy, 0.0] for j in range(ny) for i in range(nx)]
    else:
        pts = [[i*dx, j*dy]      for j in range(ny) for i in range(nx)]
    return np.array(pts, dtype=float)


def _dim_template(dim):
    
    return "Vec3d" if dim == "3d" else "Vec2d"



# =========== Q1 quad elements =========================


class _QuadElement:
    LABEL = "Q1 quad"

    @staticmethod
    def get_quads(nx, ny):
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
    def _shape_q1(xi, eta):
        N = 0.25 * np.array([
            (1 - xi)*(1 - eta),
            (1 + xi)*(1 - eta),
            (1 + xi)*(1 + eta),
            (1 - xi)*(1 + eta),
        ])
        dN_dxi  = 0.25 * np.array([-(1-eta),  (1-eta),  (1+eta), -(1+eta)])
        dN_deta = 0.25 * np.array([-(1-xi),  -(1+xi),   (1+xi),  (1-xi)])
        return N, dN_dxi, dN_deta

    @staticmethod
    def compute_nodal_forces(nodes_2d, L, E, nx, ny):
        """Always uses the first two coordinate columns (x, y)."""
        quads = _QuadElement.get_quads(nx, ny)
        # Work with xy columns regardless of 2d/3d storage
        xy = nodes_2d[:, :2]
        F  = np.zeros((len(xy), 2))

        # Body forces -- Gauss quadrature over each quad
        for quad in quads:
            xe, ye = xy[quad, 0], xy[quad, 1]
            for xi, wi in zip(GAUSS_PTS, GAUSS_WTS):
                for eta, wj in zip(GAUSS_PTS, GAUSS_WTS):
                    N, dN_dxi, dN_deta = _QuadElement._shape_q1(xi, eta)
                    J    = np.array([[dN_dxi@xe, dN_dxi@ye],
                                     [dN_deta@xe, dN_deta@ye]])
                    detJ = np.linalg.det(J)
                    xg, yg = N@xe, N@ye
                    fx, fy = fx_body(xg, yg, E, L), fy_body(xg, yg, E, L)
                    w = wi * wj * detJ
                    for a, node in enumerate(quad):
                        F[node, 0] += N[a] * fx * w
                        F[node, 1] += N[a] * fy * w

        # Neumann BC — bottom edge (y=0)
        for i in range(nx - 1):
            k0, k1 = i, i + 1
            x0, y0 = xy[k0]; x1, y1 = xy[k1]
            Le = np.hypot(x1-x0, y1-y0)
            for xi, wi in zip(GAUSS_PTS, GAUSS_WTS):
                t = 0.5*(xi + 1)
                xg, yg = (1-t)*x0 + t*x1, (1-t)*y0 + t*y1
                Tx, Ty = traction(xg, yg, 0., -1., E, L)
                N0, N1 = 0.5*(1-xi), 0.5*(1+xi)
                w = wi * Le / 2.0
                F[k0, 0] += N0*Tx*w; F[k0, 1] += N0*Ty*w
                F[k1, 0] += N1*Tx*w; F[k1, 1] += N1*Ty*w

        # Neumann BC — top edge (y=L)
        for i in range(nx - 1):
            k0 = (ny-1)*nx + i
            k1 = (ny-1)*nx + i + 1
            x0, y0 = xy[k0]; x1, y1 = xy[k1]
            Le = np.hypot(x1-x0, y1-y0)
            for xi, wi in zip(GAUSS_PTS, GAUSS_WTS):
                t = 0.5*(xi + 1)
                xg, yg = (1-t)*x0 + t*x1, (1-t)*y0 + t*y1
                Tx, Ty = traction(xg, yg, 0., +1., E, L)
                N0, N1 = 0.5*(1-xi), 0.5*(1+xi)
                w = wi * Le / 2.0
                F[k0, 0] += N0*Tx*w; F[k0, 1] += N0*Ty*w
                F[k1, 0] += N1*Tx*w; F[k1, 1] += N1*Ty*w
        return F  # shape (N, 2) — z-component added in create_sofa_scene if needed

    @staticmethod
    def create_sofa_scene(rootNode, L=1.0, E=1e6, nx=10, ny=10,
                          with_visual=True, dim="2d"):
        """
        Build SOFA scene for Q1 quads.
        dim : "2d"  →  Vec2d template
              "3d"  →  Vec3d template, z coordinate fixed at 0
        Returns (dofs, nodes_2d, quads).
        """
        tmpl = _dim_template(dim)

        rootNode.addObject("RequiredPlugin", pluginName="Sofa.Component.Visual")
        rootNode.addObject("RequiredPlugin", pluginName=[
            "Elasticity", "Sofa.Component.Constraint.Projective",
            "Sofa.Component.Engine.Select", "Sofa.Component.LinearSolver.Direct",
            "Sofa.Component.MechanicalLoad", "Sofa.Component.ODESolver.Backward",
            "Sofa.Component.StateContainer", "Sofa.Component.Topology.Container.Dynamic",
        ])
        rootNode.addObject("DefaultAnimationLoop")
        if with_visual:
            rootNode.addObject("VisualStyle",
                               displayFlags="showBehaviorModels showForceFields")

        nodes_2d = get_nodes_2d(L, nx, ny, dim=dim)
        quads_np = _QuadElement.get_quads(nx, ny)

        Beam = rootNode.addChild("Beam")
        Beam.addObject("StaticSolver", name="staticSolver", printLog=False)
        Beam.addObject('NewtonRaphsonSolver',
                  name="newtonSolver",
                  maxNbIterationsNewton=1,
                  absoluteResidualStoppingThreshold=1e-10,
                  printLog=False)
        Beam.addObject("SparseLDLSolver", name="linearSolver",
                       template="CompressedRowSparseMatrixd")

        dofs = Beam.addObject("MechanicalObject", name="dofs", template=tmpl,
                              position=nodes_2d.tolist(), showObject=with_visual,
                              showObjectScale=0.005*L)

        Beam.addObject("QuadSetTopologyContainer", name="topology",
                       quads=quads_np.tolist())
        Beam.addObject("QuadSetTopologyModifier")
        Beam.addObject("LinearSmallStrainFEMForceField", name="FEM", template=tmpl,
                       youngModulus=E, poissonRatio=0.0, topology="@topology")

        e = 1e-4 * L / max(nx-1, ny-1)
        if dim == "3d":
            box_left  = [-e, -e, -e, e,   L+e, e]
            box_right = [L-e, -e, -e, L+e, L+e, e]
        else:
            box_left  = [-e, -e, 0, e,   L+e, 0]
            box_right = [L-e, -e, 0, L+e, L+e, 0]

        Beam.addObject("BoxROI", name="box_left", template=tmpl,
                       box=box_left, drawBoxes=with_visual)
        Beam.addObject("FixedProjectiveConstraint", name="fix_left", template=tmpl,
                       indices="@box_left.indices")
        Beam.addObject("BoxROI", name="box_right", template=tmpl,
                       box=box_right, drawBoxes=with_visual)
        Beam.addObject("FixedProjectiveConstraint", name="fix_right", template=tmpl,
                       indices="@box_right.indices")

        # Build nodal force vector; pad with z=0 column for Vec3d
        F_xy = _QuadElement.compute_nodal_forces(nodes_2d, L, E, nx, ny)
        if dim == "3d":
            F_all = np.hstack([F_xy, np.zeros((len(F_xy), 1))])
        else:
            F_all = F_xy

        idx = " ".join(str(k) for k in range(len(nodes_2d)))
        if dim == "3d":
            frc = " ".join(f"{F_all[k,0]} {F_all[k,1]} {F_all[k,2]}"
                           for k in range(len(nodes_2d)))
        else:
            frc = " ".join(f"{F_all[k,0]} {F_all[k,1]}"
                           for k in range(len(nodes_2d)))

        Beam.addObject("ConstantForceField", name="MMS_forces", template=tmpl,
                       indices=idx, forces=frc)
        return dofs, nodes_2d, quads_np

    @staticmethod
    def compute_l2(nodes_2d, ux, uy, L, nx, ny, quads_np):
        xy   = nodes_2d[:, :2]
        err2 = 0.0
        for quad in quads_np:
            xe, ye = xy[quad, 0], xy[quad, 1]
            for xi, wi in zip(GAUSS_PTS, GAUSS_WTS):
                for eta, wj in zip(GAUSS_PTS, GAUSS_WTS):
                    N, dN_dxi, dN_deta = _QuadElement._shape_q1(xi, eta)
                    J    = np.array([[dN_dxi@xe, dN_dxi@ye],
                                     [dN_deta@xe, dN_deta@ye]])
                    detJ = np.linalg.det(J)
                    xg, yg = N@xe, N@ye
                    ux_h = N @ ux[quad]; uy_h = N @ uy[quad]
                    ex = ux_h - ux_mms(xg, yg, L)
                    ey = uy_h - uy_mms(xg, yg, L)
                    err2 += (ex**2 + ey**2) * wi * wj * detJ
        return np.sqrt(err2)

    @staticmethod
    def compute_h1(nodes_2d, ux, uy, L, quads_np, **_ignored):
        xy   = nodes_2d[:, :2]
        err2 = 0.0
        for quad in quads_np:
            xe, ye = xy[quad, 0], xy[quad, 1]
            for xi, wi in zip(GAUSS_PTS, GAUSS_WTS):
                for eta, wj in zip(GAUSS_PTS, GAUSS_WTS):
                    N, dN_dxi, dN_deta = _QuadElement._shape_q1(xi, eta)
                    J    = np.array([[dN_dxi@xe, dN_dxi@ye],
                                     [dN_deta@xe, dN_deta@ye]])
                    detJ = np.linalg.det(J)
                    Jinv = np.linalg.inv(J)
                    xg, yg = N@xe, N@ye
                    dN_dx = Jinv[0,0]*dN_dxi + Jinv[1,0]*dN_deta
                    dN_dy = Jinv[0,1]*dN_dxi + Jinv[1,1]*dN_deta
                    dux_dx_h = dN_dx @ ux[quad]; dux_dy_h = dN_dy @ ux[quad]
                    duy_dx_h = dN_dx @ uy[quad]; duy_dy_h = dN_dy @ uy[quad]
                    err2 += (
                        (dux_dx_h - dux_dx(xg,yg,L))**2 +
                        (dux_dy_h - dux_dy(xg,yg,L))**2 +
                        (duy_dx_h - duy_dx(xg,yg,L))**2 +
                        (duy_dy_h - duy_dy(xg,yg,L))**2
                    ) * wi * wj * detJ
        return np.sqrt(err2)


# P1 triangle elements


class _TriElement:
    LABEL = "P1 tri"

    @staticmethod
    def get_triangles(nx, ny):
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
    def compute_nodal_forces(nodes_2d, L, E, nx, ny):
        """Always uses the first two coordinate columns (x, y)."""
        dx, dy = L / (nx - 1), L / (ny - 1)
        xy  = nodes_2d[:, :2]
        F   = np.zeros((len(xy), 2))
        eps = 1e-10

        # Body forces — trapezoidal nodal quadrature
        for k, (xk, yk) in enumerate(xy):
            i = round(xk / dx); j = round(yk / dy)
            wx = dx/2 if (i == 0 or i == nx-1) else dx
            wy = dy/2 if (j == 0 or j == ny-1) else dy
            F[k, 0] += fx_body(xk, yk, E, L) * wx * wy
            F[k, 1] += fy_body(xk, yk, E, L) * wx * wy

        # Neumann BC — bottom edge (y=0)
        for i in range(nx):
            k = i
            xk, yk = xy[k]
            if xk < eps or xk > L - eps: continue
            wx = dx/2 if (i == 0 or i == nx-1) else dx
            Tx, Ty = traction(xk, yk, 0., -1., E, L)
            F[k, 0] += Tx * wx; F[k, 1] += Ty * wx

        # Neumann BC — top edge (y=L)
        for i in range(nx):
            k = (ny-1)*nx + i
            xk, yk = xy[k]
            if xk < eps or xk > L - eps: continue
            wx = dx/2 if (i == 0 or i == nx-1) else dx
            Tx, Ty = traction(xk, yk, 0., +1., E, L)
            F[k, 0] += Tx * wx; F[k, 1] += Ty * wx
        return F  # shape (N, 2)

    @staticmethod
    def create_sofa_scene(rootNode, L=1.0, E=1e6, nx=10, ny=10,
                          with_visual=True, dim="2d"):
     
        tmpl = _dim_template(dim)

        rootNode.addObject("RequiredPlugin", pluginName="Sofa.Component.Visual")
        rootNode.addObject("RequiredPlugin", pluginName=[
            "Elasticity", "Sofa.Component.Constraint.Projective",
            "Sofa.Component.Engine.Select", "Sofa.Component.LinearSolver.Direct",
            "Sofa.Component.MechanicalLoad", "Sofa.Component.ODESolver.Backward",
            "Sofa.Component.StateContainer", "Sofa.Component.Topology.Container.Dynamic",
        ])
        rootNode.addObject("DefaultAnimationLoop")
        if with_visual:
            rootNode.addObject("VisualStyle",
                               displayFlags="showBehaviorModels showForceFields")

        nodes_2d = get_nodes_2d(L, nx, ny, dim=dim)
        tris_np  = _TriElement.get_triangles(nx, ny)

        Beam = rootNode.addChild("Beam")
        Beam.addObject("StaticSolver", name="staticSolver", printLog=False)
        Beam.addObject('NewtonRaphsonSolver',
                  name="newtonSolver",
                  maxNbIterationsNewton=1,
                  absoluteResidualStoppingThreshold=1e-10,
                  printLog=False)
        Beam.addObject("SparseLDLSolver", name="linearSolver",
                       template="CompressedRowSparseMatrixd")

        dofs = Beam.addObject("MechanicalObject", name="dofs", template=tmpl,
                              position=nodes_2d.tolist(), showObject=with_visual,
                              showObjectScale=0.005*L)

        Beam.addObject("TriangleSetTopologyContainer", name="topology",
                       triangles=tris_np.tolist())
        Beam.addObject("TriangleSetTopologyModifier")
        Beam.addObject("LinearSmallStrainFEMForceField", name="FEM", template=tmpl,
                       youngModulus=E, poissonRatio=0.0, topology="@topology")

        e = 1e-4 * L / max(nx-1, ny-1)
        if dim == "3d":
            box_left  = [-e, -e, -e, e,   L+e, e]
            box_right = [L-e, -e, -e, L+e, L+e, e]
        else:
            box_left  = [-e, -e, 0, e,   L+e, 0]
            box_right = [L-e, -e, 0, L+e, L+e, 0]

        Beam.addObject("BoxROI", name="box_left", template=tmpl,
                       box=box_left, drawBoxes=with_visual)
        Beam.addObject("FixedProjectiveConstraint", name="fix_left", template=tmpl,
                       indices="@box_left.indices")
        Beam.addObject("BoxROI", name="box_right", template=tmpl,
                       box=box_right, drawBoxes=with_visual)
        Beam.addObject("FixedProjectiveConstraint", name="fix_right", template=tmpl,
                       indices="@box_right.indices")

        # Build nodal force vector; pad with z=0 column for Vec3d
        F_xy = _TriElement.compute_nodal_forces(nodes_2d, L, E, nx, ny)
        if dim == "3d":
            F_all = np.hstack([F_xy, np.zeros((len(F_xy), 1))])
        else:
            F_all = F_xy

        idx = " ".join(str(k) for k in range(len(nodes_2d)))
        if dim == "3d":
            frc = " ".join(f"{F_all[k,0]} {F_all[k,1]} {F_all[k,2]}"
                           for k in range(len(nodes_2d)))
        else:
            frc = " ".join(f"{F_all[k,0]} {F_all[k,1]}"
                           for k in range(len(nodes_2d)))

        Beam.addObject("ConstantForceField", name="MMS_forces", template=tmpl,
                       indices=idx, forces=frc)
        return dofs, nodes_2d, tris_np

    @staticmethod
    def compute_l2(nodes_2d, ux, uy, L, nx, ny, tris_np):
        """L2 error using centroid quadrature."""
        dx, dy = L / (nx - 1), L / (ny - 1)
        area   = dx * dy / 2.0
        xy     = nodes_2d[:, :2]
        err2   = 0.0
        for i0, i1, i2 in tris_np:
            xc = (xy[i0,0] + xy[i1,0] + xy[i2,0]) / 3
            yc = (xy[i0,1] + xy[i1,1] + xy[i2,1]) / 3
            ux_c = (ux[i0] + ux[i1] + ux[i2]) / 3
            uy_c = (uy[i0] + uy[i1] + uy[i2]) / 3
            err2 += ((ux_c - ux_mms(xc,yc,L))**2 +
                     (uy_c - uy_mms(xc,yc,L))**2) * area
        return np.sqrt(err2)

    @staticmethod
    def _grad_p1(nodes_2d, u, tri):
        xy = nodes_2d[:, :2]
        i0, i1, i2 = tri
        x0, y0 = xy[i0]; x1, y1 = xy[i1]; x2, y2 = xy[i2]
        A2 = (x1-x0)*(y2-y0) - (x2-x0)*(y1-y0)
        dphi0_dx = (y1-y2)/A2; dphi0_dy = (x2-x1)/A2
        dphi1_dx = (y2-y0)/A2; dphi1_dy = (x0-x2)/A2
        dphi2_dx = (y0-y1)/A2; dphi2_dy = (x1-x0)/A2
        dudx = u[i0]*dphi0_dx + u[i1]*dphi1_dx + u[i2]*dphi2_dx
        dudy = u[i0]*dphi0_dy + u[i1]*dphi1_dy + u[i2]*dphi2_dy
        return dudx, dudy, abs(A2)/2.0

    @staticmethod
    def compute_h1(nodes_2d, ux, uy, L, tris_np, **_ignored):
        err2 = 0.0
        for tri in tris_np:
            i0, i1, i2 = tri
            xy = nodes_2d[:, :2]
            xc = (xy[i0,0] + xy[i1,0] + xy[i2,0]) / 3
            yc = (xy[i0,1] + xy[i1,1] + xy[i2,1]) / 3
            dux_dx_h, dux_dy_h, area = _TriElement._grad_p1(nodes_2d, ux, tri)
            duy_dx_h, duy_dy_h, _    = _TriElement._grad_p1(nodes_2d, uy, tri)
            err2 += (
                (dux_dx_h - dux_dx(xc,yc,L))**2 +
                (dux_dy_h - dux_dy(xc,yc,L))**2 +
                (duy_dx_h - duy_dx(xc,yc,L))**2 +
                (duy_dy_h - duy_dy(xc,yc,L))**2
            ) * area
        return np.sqrt(err2)


# ─────────────────────────────────────────────────────────────────────────────
# Simulation runner
# ─────────────────────────────────────────────────────────────────────────────

def run_simulation(elem, L, E, nx, ny, dim="2d"):
    root = Sofa.Core.Node("root")
    dofs, nodes_2d, conn = elem.create_sofa_scene(
        root, L=L, E=E, nx=nx, ny=ny, with_visual=False, dim=dim
    )
    Sofa.Simulation.init(root)
    pos0 = dofs.position.array().copy()
    Sofa.Simulation.animate(root, root.dt.value)
    pos1 = dofs.position.array().copy()
    Sofa.Simulation.unload(root)
    # Extract x and y displacement columns (works for both Vec2d and Vec3d)
    ux = pos1[:, 0] - pos0[:, 0]
    uy = pos1[:, 1] - pos0[:, 1]
    return nodes_2d, conn, ux, uy


# ─────────────────────────────────────────────────────────────────────────────
# Convergence study
# ─────────────────────────────────────────────────────────────────────────────

def convergence_study(elem_specs, L, E, nx_values, dim="2d", results_dir=None):
    """
    Run convergence study for each element type in elem_specs.
    elem_specs : list of dicts with keys 'elem', 'label', 'marker', 'color'
    dim        : "2d" or "3d"
    """
    if results_dir is None:
        results_dir = RESULTS_DIR
    os.makedirs(results_dir, exist_ok=True)

    fig_l2, ax_l2 = plt.subplots(figsize=(7, 5))
    fig_h1, ax_h1 = plt.subplots(figsize=(7, 5))

    for spec in elem_specs:
        elem, label = spec["elem"], spec["label"]
        marker, color = spec.get("marker", "o"), spec.get("color", None)

        hs, errs_l2, errs_h1 = [], [], []
        hdr = (f"{'nx':>5} | {'h':>8} | {'L2':>14} | {'ord_L2':>7} "
               f"| {'H1':>14} | {'ord_H1':>7}")
        tag      = label.replace(" ", "_")
        txt_path = os.path.join(results_dir, f"convergence_{tag}_{dim}.txt")

        print(f"\n── Convergence  {label}  [{dim.upper()}] ──\n{hdr}")
        with open(txt_path, "w") as f:
            f.write(f"Convergence MMS 2D — {label} [{dim.upper()}]\n{hdr}\n")
            for k, nx in enumerate(nx_values):
                ny = nx
                h  = L / (nx - 1)
                nodes_2d, conn, ux, uy = run_simulation(elem, L, E, nx, ny, dim=dim)
                l2 = elem.compute_l2(nodes_2d, ux, uy, L, nx, ny, conn)
                h1 = elem.compute_h1(nodes_2d, ux, uy, L, conn)
                hs.append(h); errs_l2.append(l2); errs_h1.append(h1)

                ord_l2 = (f"{np.log(l2/errs_l2[k-1])/np.log(h/hs[k-1]):.2f}"
                          if k > 0 else "")
                ord_h1 = (f"{np.log(h1/errs_h1[k-1])/np.log(h/hs[k-1]):.2f}"
                          if k > 0 else "")
                line = (f"{nx:5d} | {h:8.4f} | {l2:14.6e} | {ord_l2:>7} "
                        f"| {h1:14.6e} | {ord_h1:>7}")
                print(line); f.write(line + "\n")

        hs_a, l2_a, h1_a = np.array(hs), np.array(errs_l2), np.array(errs_h1)
        kw = dict(lw=2, ms=7, color=color)
        ax_l2.loglog(hs_a, l2_a, f"{marker}-",  label=f"{label} [{dim}]", **kw)
        ax_h1.loglog(hs_a, h1_a, f"{marker}--", label=f"{label} [{dim}]", **kw)

    # Reference slopes
    h_ref = np.array([hs_a[0], hs_a[-1]])
    ax_l2.loglog(h_ref, l2_a[0]*(h_ref/hs_a[0])**2, ":", color="gray",
                 lw=1.2, label="O(h²)")
    ax_h1.loglog(h_ref, h1_a[0]*(h_ref/hs_a[0])**1, ":", color="gray",
                 lw=1.2, label="O(h¹)")

    for ax, ylabel, title, fname in [
        (ax_l2, "L² error",    f"L² convergence MMS 2D [{dim}]",
         f"convergence_L2_{dim}.png"),
        (ax_h1, "H¹ semi-norm", f"H¹ convergence MMS 2D [{dim}]",
         f"convergence_H1_{dim}.png"),
    ]:
        ax.set_xlabel("h"); ax.set_ylabel(ylabel); ax.set_title(title)
        ax.legend(); ax.grid(True, alpha=0.3, which="both")

    fig_l2.tight_layout(); fig_h1.tight_layout()
    fig_l2.savefig(os.path.join(results_dir, f"convergence_L2_{dim}.png"), dpi=150)
    fig_h1.savefig(os.path.join(results_dir, f"convergence_H1_{dim}.png"), dpi=150)
    plt.close(fig_l2); plt.close(fig_h1)


# ─────────────────────────────────────────────────────────────────────────────
# Point simulation (solution plots)
# ─────────────────────────────────────────────────────────────────────────────

def simulation_ponctuelle(elem, L, E, nx, ny, dim="2d", results_dir=None):
    if results_dir is None:
        results_dir = RESULTS_DIR
    os.makedirs(results_dir, exist_ok=True)

    nodes_2d, conn, ux, uy = run_simulation(elem, L, E, nx, ny, dim=dim)
    l2 = elem.compute_l2(nodes_2d, ux, uy, L, nx, ny, conn)
    h1 = elem.compute_h1(nodes_2d, ux, uy, L, conn)

    xy     = nodes_2d[:, :2]
    ux_ref = ux_mms(xy[:, 0], xy[:, 1], L)
    uy_ref = uy_mms(xy[:, 0], xy[:, 1], L)

    label = elem.LABEL
    print(f"\n[{label}]  dim={dim}  nx={nx} ny={ny}  L={L}")
    print(f"  L2             = {l2:.4e}")
    print(f"  H1 semi-norm   = {h1:.4e}")
    print(f"  max|ux-ux_mms| = {np.max(np.abs(ux - ux_ref)):.4e}")
    print(f"  max|uy-uy_mms| = {np.max(np.abs(uy - uy_ref)):.4e}")

    # 1-D profile along mid-height line
    mid_j = (ny - 1) // 2
    sl    = slice(mid_j * nx, mid_j * nx + nx)
    yc    = xy[mid_j * nx, 1]
    xf    = np.linspace(0, L, 300)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    for ax, u_sofa, u_fn, lbl, fmt in zip(
        axes,
        [ux[sl], uy[sl]],
        [lambda x: ux_mms(x, yc, L), lambda x: uy_mms(x, yc, L)],
        [r"$u_x$", r"$u_y$"],
        ["o-", "s-"],
    ):
        ax.plot(xy[sl, 0], u_sofa, fmt, color="tab:green",
                label=f"SOFA {label} [{dim}]", ms=5)
        ax.plot(xf, u_fn(xf), "--", color="tab:blue", label="MMS exact")
        ax.set_xlabel("x"); ax.set_ylabel(lbl)
        ax.legend(); ax.grid(True, alpha=0.3)
    plt.suptitle(f"MMS 2D — {label} [{dim}]  nx={nx}  |L2={l2:.2e}  H1={h1:.2e}")
    plt.tight_layout()
    tag = label.replace(" ", "_")
    plt.savefig(os.path.join(results_dir,
                             f"solution_{tag}_{dim}_nx{nx}.png"), dpi=150)
    plt.close()

    # 2-D colour maps
    tris_plot = elem.to_triangles(conn)
    x, y = xy[:, 0], xy[:, 1]
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    for ax, data, title, cmap in [
        (axes[0,0], ux,                r"$u_x$ SOFA",               "gray"),
        (axes[0,1], ux_ref,            r"$u_x$ MMS",                "gray"),
        (axes[0,2], np.abs(ux-ux_ref), r"$|u_x - u_x^{MMS}|$",     "gray_r"),
        (axes[1,0], uy,                r"$u_y$ SOFA",               "gray"),
        (axes[1,1], uy_ref,            r"$u_y$ MMS",                "gray"),
        (axes[1,2], np.abs(uy-uy_ref), r"$|u_y - u_y^{MMS}|$",     "gray_r"),
    ]:
        tc = ax.tricontourf(x, y, tris_plot.tolist(), data, levels=20, cmap=cmap)
        ax.triplot(x, y, tris_plot.tolist(), "k-", lw=0.3, alpha=0.4)
        plt.colorbar(tc, ax=ax, shrink=0.8)
        ax.set_title(title); ax.set_aspect("equal")
        ax.set_xlabel("x"); ax.set_ylabel("y")
    plt.suptitle(f"Champs 2D — {label} [{dim}]  nx={nx}")
    plt.tight_layout()
    plt.savefig(os.path.join(results_dir,
                             f"champs2D_{tag}_{dim}_nx{nx}.png"), dpi=150)
    plt.close()


# ─────────────────────────────────────────────────────────────────────────────
# Element instances
# ─────────────────────────────────────────────────────────────────────────────

element_quad = _QuadElement()
element_tri  = _TriElement()


if __name__ == "__main__":

    
    #   DIM = "2d"     Vec2d template (plan)  
    #   DIM = "3d"     Vec3d template (z=0)   

    DIM = "2d"       

    L, E    = 1.0, 1e6

    nx_vals = [10,20,30,40,50,60]



    specs = [
        {"elem": element_quad, "label": "Q1 quad", "marker": "o", "color": "C0"},
        {"elem": element_tri,  "label": "P1 tri",  "marker": "s", "color": "C1"},
    ]

    convergence_study(specs, L, E, nx_vals, dim=DIM)
    simulation_ponctuelle(element_quad, L, E, nx=20, ny=20, dim=DIM)
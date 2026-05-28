import json
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

import Sofa
import Sofa.Core
import Sofa.Simulation
import SofaRuntime

# Make the parent MMS/ directory importable so we can pull in fem.py.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from fem import hex_q1_rule              # re-exported for case files
from elements import element_hex          # re-exported for case files
from solid_solution import SolidSolution3D
from output import write_solution_table
from scene import NodalForceAssembler

RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")


# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

def load_params(path=None):
    if path is None:
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "params.json")
    with open(path) as f:
        return json.load(f)


# ============== Mesh helpers ============================

def get_nodes_3d(L, nx, ny, nz):
    dx = L / (nx - 1)
    dy = L / (ny - 1)
    dz = L / (nz - 1)
    pts = [[i*dx, j*dy, k*dz]
           for k in range(nz) for j in range(ny) for i in range(nx)]
    return np.array(pts, dtype=float)


# ---------------------------------------------------------------------------
# SOFA scene
# ---------------------------------------------------------------------------

def _solid_force_compute(element, mms, L, E, nu, nx, ny, nz):
    """Return a `(nodes, topology) -> (N, 3) force array` for NodalForceAssembler."""
    def compute(nodes, topology):
        conn = element.read_connectivity(topology)
        return element.compute_nodal_forces(
            nodes, conn, mms, L, E, nu, nx, ny, nz)
    return compute

def build_solid_scene(rootNode, mms, element, L=1.0, E=1e6, nu=0.3,
                      nx=6, ny=6, nz=6, with_visual=True):
    """Build a SOFA scene for `mms` on the 3D `element` strategy.

    Returns (dofs, topology). Nodes and connectivity become available
    after `Sofa.Simulation.init(root)` runs, via
    `dofs.rest_position.array()` and `element.read_connectivity(topology)`.
    """
    rootNode.addObject("RequiredPlugin", pluginName=[
        "Elasticity",
        "Sofa.Component.Constraint.Projective",
        "Sofa.Component.Engine.Select",
        "Sofa.Component.LinearSolver.Direct",
        "Sofa.Component.MechanicalLoad",
        "Sofa.Component.ODESolver.Backward",
        "Sofa.Component.StateContainer",
        "Sofa.Component.Topology.Container.Grid",
        "Sofa.Component.Topology.Container.Dynamic",
        "Sofa.Component.Topology.Mapping",
        "Sofa.Component.Visual",
    ])
    rootNode.addObject("DefaultAnimationLoop")
    if with_visual:
        rootNode.addObject("VisualStyle",
                           displayFlags="showBehaviorModels showForceFields")

    nodes_3d = get_nodes_3d(L, nx, ny, nz)

    # Grid must be added *before* Solid: SOFA inits children in insertion
    # order, and Solid's topology container resolves `@../Grid/grid.position`
    # during its own init — so Grid has to be initialised first.
    Grid = rootNode.addChild("Grid")
    Grid.addObject("RegularGridTopology", name="grid",
                   nx=nx, ny=ny, nz=nz,
                   min=[0.0, 0.0, 0.0], max=[L, L, L])

    Solid = rootNode.addChild("Solid")
    Solid.addObject("StaticSolver", name="staticSolver", printLog=False)
    Solid.addObject("NewtonRaphsonSolver", name="newtonSolver",
                    maxNbIterationsNewton=1,
                    absoluteResidualStoppingThreshold=1e-10,
                    printLog=False)
    Solid.addObject("SparseLDLSolver", name="linearSolver",
                    template="CompressedRowSparseMatrixd")

    dofs = Solid.addObject("MechanicalObject", name="dofs", template="Vec3d",
                           position=nodes_3d.tolist(),
                           showObject=with_visual, showObjectScale=0.005 * L)

    topology = element.add_topology(Solid)

    Solid.addObject("LinearSmallStrainFEMForceField", name="FEM", template="Vec3d",
                    youngModulus=E, poissonRatio=nu, topology="@topology")

    mms.apply_bcs(Solid, nodes_3d, L)

    # Placeholder force field filled in by the controller after init.
    n_nodes = len(nodes_3d)
    force_field = Solid.addObject("ConstantForceField", name="MMS_forces",
                                  template="Vec3d",
                                  indices=list(range(n_nodes)),
                                  forces=[[0.0, 0.0, 0.0]] * n_nodes)

    Solid.addObject(NodalForceAssembler(
        dofs=dofs, topology=topology, force_field=force_field,
        compute_forces=_solid_force_compute(element, mms, L, E, nu, nx, ny, nz),
        name="nodalForceAssembler"))

    return dofs, topology


# ─────────────────────────────────────────────────────────────────────────────
# Simulation runner
# ─────────────────────────────────────────────────────────────────────────────

def solve_solid(elem, mms, L, E, nu, nx, ny, nz):
    """Build, init, and run one static step. Returns a SolidSolution3D snapshot."""
    root = Sofa.Core.Node("root")
    dofs, topology = build_solid_scene(
        root, mms, elem, L=L, E=E, nu=nu,
        nx=nx, ny=ny, nz=nz, with_visual=False
    )
    Sofa.Simulation.init(root)
    nodes_3d = dofs.rest_position.array().copy()
    conn     = elem.read_connectivity(topology)
    pos0     = dofs.position.array().copy()
    Sofa.Simulation.animate(root, root.dt.value)
    pos1     = dofs.position.array().copy()
    Sofa.Simulation.unload(root)
    ux = pos1[:, 0] - pos0[:, 0]
    uy = pos1[:, 1] - pos0[:, 1]
    uz = pos1[:, 2] - pos0[:, 2]
    return SolidSolution3D(nodes=nodes_3d, conn=conn, ux=ux, uy=uy, uz=uz)


# ---------------------------------------------------------------------------
# Output helpers (mirror 2D beam.py)
# ---------------------------------------------------------------------------

def plot_solution_profile(stem, sol, mms, L, nx, ny, nz, label, nu, l2, h1):
    """Save 1-D centerline profiles (ux(x), uy(y), uz(z)) to results/<stem>.png."""
    os.makedirs(RESULTS_DIR, exist_ok=True)
    xyz = sol.nodes[:, :3]

    mid_i, mid_j, mid_k = nx // 2, ny // 2, nz // 2

    def nidx(i, j, k):
        return i + j * nx + k * nx * ny

    x_fine = np.linspace(0, L, 200)

    line_x = [nidx(i, mid_j, mid_k) for i in range(nx)]
    line_y = [nidx(mid_i, j, mid_k) for j in range(ny)]
    line_z = [nidx(mid_i, mid_j, k) for k in range(nz)]

    yc, zc = xyz[line_x[0], 1], xyz[line_x[0], 2]
    xc2, zc2 = xyz[line_y[0], 0], xyz[line_y[0], 2]
    xc3, yc3 = xyz[line_z[0], 0], xyz[line_z[0], 1]

    ux_ex, _, _ = mms.u_ex(x_fine,
                           np.full_like(x_fine, yc),
                           np.full_like(x_fine, zc), L)
    _, uy_ex, _ = mms.u_ex(np.full_like(x_fine, xc2),
                           x_fine,
                           np.full_like(x_fine, zc2), L)
    _, _, uz_ex = mms.u_ex(np.full_like(x_fine, xc3),
                           np.full_like(x_fine, yc3),
                           x_fine, L)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for ax, coord, sofa, exact, ylabel, axname in [
        (axes[0], xyz[line_x, 0], sol.ux[line_x], ux_ex, r"$u_x$", "x"),
        (axes[1], xyz[line_y, 1], sol.uy[line_y], uy_ex, r"$u_y$", "y"),
        (axes[2], xyz[line_z, 2], sol.uz[line_z], uz_ex, r"$u_z$", "z"),
    ]:
        ax.plot(coord,  sofa,  "o-", color="tab:green",
                label=f"SOFA {label}", ms=5)
        ax.plot(x_fine, exact, "--", color="tab:blue", label="MMS exact")
        ax.set_xlabel(axname); ax.set_ylabel(ylabel)
        ax.legend(); ax.grid(True, alpha=0.3)
    fig.suptitle(f"{mms.name} — {label}  "
                 f"nu={nu}  nx={nx}  |L2={l2:.2e}  H1={h1:.2e}")
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, f"{stem}.png"), dpi=150)
    plt.close(fig)


def plot_solution_slices(stem, sol, mms, L, nx, ny, nz, label, nu):
    """Save midplane heatmaps (xy at z=mid, xz at y=mid, yz at x=mid) for ux/uy/uz.

    Mirrors the 2-row × 3-column layout of the existing 3D plot_displacement:
    each column is a slice plane, each row a displacement component (two per plane).
    """
    os.makedirs(RESULTS_DIR, exist_ok=True)
    xyz = sol.nodes[:, :3]
    mid_i, mid_j, mid_k = nx // 2, ny // 2, nz // 2

    def nidx(i, j, k):
        return i + j * nx + k * nx * ny

    sl_z = np.array([[nidx(i, j, mid_k) for i in range(nx)] for j in range(ny)])
    sl_y = np.array([[nidx(i, mid_j, k) for i in range(nx)] for k in range(nz)])
    sl_x = np.array([[nidx(mid_i, j, k) for j in range(ny)] for k in range(nz)])

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    # Column 0: z = mid plane — show ux, uy
    X, Y = xyz[sl_z, 0], xyz[sl_z, 1]
    for row, u, title in [(0, sol.ux, r"$u_x$ z=mid"),
                          (1, sol.uy, r"$u_y$ z=mid")]:
        im = axes[row, 0].pcolormesh(X, Y, u[sl_z], cmap="RdBu_r", shading="auto")
        axes[row, 0].set(xlabel="x", ylabel="y", title=title)
        plt.colorbar(im, ax=axes[row, 0])

    # Column 1: y = mid plane — show ux, uz
    X, Z = xyz[sl_y, 0], xyz[sl_y, 2]
    for row, u, title in [(0, sol.ux, r"$u_x$ y=mid"),
                          (1, sol.uz, r"$u_z$ y=mid")]:
        im = axes[row, 1].pcolormesh(X, Z, u[sl_y], cmap="RdBu_r", shading="auto")
        axes[row, 1].set(xlabel="x", ylabel="z", title=title)
        plt.colorbar(im, ax=axes[row, 1])

    # Column 2: x = mid plane — show uy, uz
    Y, Z = xyz[sl_x, 1], xyz[sl_x, 2]
    for row, u, title in [(0, sol.uy, r"$u_y$ x=mid"),
                          (1, sol.uz, r"$u_z$ x=mid")]:
        im = axes[row, 2].pcolormesh(Y, Z, u[sl_x], cmap="RdBu_r", shading="auto")
        axes[row, 2].set(xlabel="y", ylabel="z", title=title)
        plt.colorbar(im, ax=axes[row, 2])

    fig.suptitle(f"Fields 3D — {label}  {mms.name}  nu={nu}  nx={nx}")
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, f"{stem}.png"), dpi=150)
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# Single-case driver
# ─────────────────────────────────────────────────────────────────────────────

def run_reference_scene(elem, mms):
    """Solve one MMS case at the reference mesh, write the solution table and plots.

    All parameters come from params.json (top-level + `reference` block).
    """
    cfg = load_params()
    ref = cfg["reference"]
    L, E = cfg["length"], cfg["youngModulus"]
    nu   = ref["nu"]
    nx = ny = nz = ref["nx"]

    sol = solve_solid(elem, mms, L, E, nu, nx, ny, nz)
    l2  = elem.compute_l2(sol, mms, L)
    h1  = elem.compute_h1(sol, mms, L)

    label = elem.LABEL
    tag   = label.replace(" ", "_")
    stem  = f"{mms.name}_{tag}_nu{nu}_nx{nx}"

    xyz = sol.nodes[:, :3]
    write_solution_table(f"solution_{stem}", xyz,
                         np.column_stack([sol.ux, sol.uy, sol.uz]),
                         lambda xi, yi, zi: mms.u_ex(xi, yi, zi, L),
                         RESULTS_DIR, {"L2": l2, "H1_semi": h1})
    plot_solution_profile(f"solution_{stem}", sol, mms, L, nx, ny, nz,
                          label, nu, l2, h1)
    plot_solution_slices (f"fields3D_{stem}", sol, mms, L, nx, ny, nz,
                          label, nu)


def case_scene(mms, element):
    """Return a `createScene(rootNode)` bound to this MMS and element type.

    All parameters come from params.json (top-level + `reference` block). Each
    case file exposes:
        createScene = case_scene(mms, element_hex)
    so that `runSofa sinus_neumann.py` loads the default scene.
    """
    def createScene(rootNode):
        cfg = load_params()
        ref = cfg["reference"]
        build_solid_scene(rootNode, mms, element,
                          L=cfg["length"], E=cfg["youngModulus"],
                          nu=ref["nu"],
                          nx=ref["nx"], ny=ref["nx"], nz=ref["nx"],
                          with_visual=True)
        return rootNode
    return createScene

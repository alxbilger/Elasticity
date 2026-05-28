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
from fem import (
    assemble_nodal_forces_3d,
    assemble_traction_3d,
    l2_error_3d,
    h1_semi_error_3d,
    hex_q1_rule,
)
from solid_solution import SolidSolution3D

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


def _boundary_quads(nx, ny, nz):
    """Return (xm, xp, ym, yp, zm, zp) quad lists for a structured nx×ny×nz grid.

    Each quad is a 4-tuple of node indices. Node index convention:
    `idx(i, j, k) = i + j*nx + k*nx*ny` (matches SOFA RegularGridTopology).
    Orientation is consistent within each face but its sign does not affect
    `assemble_traction_3d`, which uses |t_xi × t_eta| for the surface area
    and takes the outward normal from the caller.
    """
    def idx(i, j, k):
        return i + j * nx + k * nx * ny

    xm = [(idx(0, j, k),    idx(0, j+1, k),    idx(0, j+1, k+1),    idx(0, j, k+1))
          for k in range(nz - 1) for j in range(ny - 1)]
    xp = [(idx(nx-1, j, k), idx(nx-1, j+1, k), idx(nx-1, j+1, k+1), idx(nx-1, j, k+1))
          for k in range(nz - 1) for j in range(ny - 1)]
    ym = [(idx(i, 0, k),    idx(i+1, 0, k),    idx(i+1, 0, k+1),    idx(i, 0, k+1))
          for k in range(nz - 1) for i in range(nx - 1)]
    yp = [(idx(i, ny-1, k), idx(i+1, ny-1, k), idx(i+1, ny-1, k+1), idx(i, ny-1, k+1))
          for k in range(nz - 1) for i in range(nx - 1)]
    zm = [(idx(i, j, 0),    idx(i+1, j, 0),    idx(i+1, j+1, 0),    idx(i, j+1, 0))
          for j in range(ny - 1) for i in range(nx - 1)]
    zp = [(idx(i, j, nz-1), idx(i+1, j, nz-1), idx(i+1, j+1, nz-1), idx(i, j+1, nz-1))
          for j in range(ny - 1) for i in range(nx - 1)]
    return xm, xp, ym, yp, zm, zp


# ---------------------------------------------------------------------------
# Element strategies
#
# Each element type is a thin class declaring:
#   LABEL              : human-readable identifier (used in plot legends)
#   ELEMENT_RULE       : an element rule (e.g. hex_q1_rule(2))
#   add_topology(Solid) : wire the SOFA topology and return the topology
#       container that holds the connectivity.
#       Uses RegularGridTopology under a sibling `Grid` child.
#   read_connectivity(topology) : extract the connectivity array from the
#       SOFA topology container, post-init.
#
# Assembly / error-norm logic lives once on _ElementBase. Connectivity is
# supplied by the caller (the NodalForceAssembler controller reads it from
# SOFA after init); the element class never holds its own copy.
# ---------------------------------------------------------------------------

class _ElementBase:
    @classmethod
    def compute_nodal_forces(cls, nodes_3d, conn, mms, L, E, nu, nx, ny, nz):
        xyz = nodes_3d[:, :3]

        F = assemble_nodal_forces_3d(
            lambda x, y, z: mms.source(x, y, z, E, nu, L),
            xyz, conn, cls._source_rule(mms))

        xm, xp, ym, yp, zm, zp = _boundary_quads(nx, ny, nz)
        sides = [(xm, -1.0, 0.0, 0.0),
                 (xp, +1.0, 0.0, 0.0),
                 (ym,  0.0, -1.0, 0.0),
                 (yp,  0.0, +1.0, 0.0),
                 (zm,  0.0, 0.0, -1.0),
                 (zp,  0.0, 0.0, +1.0)]
        for quads, nrm_x, nrm_y, nrm_z in sides:
            F += assemble_traction_3d(
                lambda x, y, z, nx=nrm_x, ny=nrm_y, nz=nrm_z:
                    mms.traction(x, y, z, nx, ny, nz, E, nu, L),
                xyz, quads)
        return F

    @classmethod
    def compute_l2(cls, sol, mms, L):
        return l2_error_3d(
            sol.nodes, sol.conn,
            np.column_stack([sol.ux, sol.uy, sol.uz]),
            lambda x, y, z: mms.u_ex(x, y, z, L),
            cls.ELEMENT_RULE)

    @classmethod
    def compute_h1(cls, sol, mms, L):
        return h1_semi_error_3d(
            sol.nodes, sol.conn,
            np.column_stack([sol.ux, sol.uy, sol.uz]),
            lambda x, y, z: mms.grad_u_ex(x, y, z, L),
            cls.ELEMENT_RULE)


class _HexElement(_ElementBase):
    LABEL        = "Q1 hex"
    ELEMENT_RULE = staticmethod(hex_q1_rule(2))   # used for L²/H¹ error norms

    @staticmethod
    def _source_rule(mms):
        rule = mms.source_quadrature_hex
        if rule is None:
            raise ValueError(
                f"{type(mms).__name__}.source_quadrature_hex must be set")
        return rule

    @staticmethod
    def add_topology(Solid):
        topology = Solid.addObject("HexahedronSetTopologyContainer",
                                   name="topology",
                                   hexahedra="@../Grid/grid.hexahedra",
                                   position="@../Grid/grid.position")
        Solid.addObject("HexahedronSetTopologyModifier")
        return topology

    @staticmethod
    def read_connectivity(topology):
        return topology.hexahedra.array().copy()


# ---------------------------------------------------------------------------
# Force assembly controller
#
# Same pattern as 2D: defer assembly to onSimulationInitDoneEvent so that
# SOFA topology components (RegularGridTopology,
# HexahedronSetTopologyContainer, ...) have been initialised and connectivity
# can be read off the topology container.
# ---------------------------------------------------------------------------

class NodalForceAssembler(Sofa.Core.Controller):
    def __init__(self, dofs, topology, force_field, element, mms,
                 L, E, nu, nx, ny, nz, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.dofs        = dofs
        self.topology    = topology
        self.force_field = force_field
        self.element     = element
        self.mms         = mms
        self.L, self.E, self.nu     = L, E, nu
        self.nx, self.ny, self.nz   = nx, ny, nz

    def onSimulationInitDoneEvent(self, event):
        nodes = self.dofs.rest_position.array().copy()
        conn  = self.element.read_connectivity(self.topology)
        F     = self.element.compute_nodal_forces(
            nodes, conn, self.mms,
            self.L, self.E, self.nu, self.nx, self.ny, self.nz)
        with self.force_field.forces.writeableArray() as forces:
            forces[:] = F


# ---------------------------------------------------------------------------
# SOFA scene
# ---------------------------------------------------------------------------

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
        element=element, mms=mms,
        L=L, E=E, nu=nu, nx=nx, ny=ny, nz=nz,
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

def write_solution_table(stem, x, y, z, ux_h, uy_h, uz_h, u_ex, error_dict):
    """
    Write per-node solution table and error summary to results/<stem>.txt.

    u_ex       : callable (x, y, z) -> (ux_ex, uy_ex, uz_ex)
    error_dict : ordered dict of {label: value} for error summary lines
    """
    os.makedirs(RESULTS_DIR, exist_ok=True)
    path = os.path.join(RESULTS_DIR, f"{stem}.txt")
    with open(path, "w") as f:
        f.write(f"{'x':>10} | {'y':>10} | {'z':>10} | "
                f"{'ux_h':>15} | {'uy_h':>15} | {'uz_h':>15} | "
                f"{'ux_ex':>15} | {'uy_ex':>15} | {'uz_ex':>15} | "
                f"{'err_x':>15} | {'err_y':>15} | {'err_z':>15}\n")
        f.write("-" * 184 + "\n")
        for xi, yi, zi, uxi, uyi, uzi in zip(x, y, z, ux_h, uy_h, uz_h):
            uxe, uye, uze = u_ex(xi, yi, zi)
            f.write(f"{xi:10.4f} | {yi:10.4f} | {zi:10.4f} | "
                    f"{uxi:15.6e} | {uyi:15.6e} | {uzi:15.6e} | "
                    f"{uxe:15.6e} | {uye:15.6e} | {uze:15.6e} | "
                    f"{abs(uxi - uxe):15.6e} | {abs(uyi - uye):15.6e} | "
                    f"{abs(uzi - uze):15.6e}\n")
        f.write("\n")
        for label, val in error_dict.items():
            f.write(f"{label:12s} = {val:.6e}\n")


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
    write_solution_table(f"solution_{stem}",
                         xyz[:, 0], xyz[:, 1], xyz[:, 2],
                         sol.ux, sol.uy, sol.uz,
                         lambda xi, yi, zi: mms.u_ex(xi, yi, zi, L),
                         {"L2": l2, "H1_semi": h1})
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


# ─────────────────────────────────────────────────────────────────────────────
# Element instances
# ─────────────────────────────────────────────────────────────────────────────

element_hex = _HexElement()

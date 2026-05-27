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
    assemble_nodal_forces_2d,
    assemble_traction_2d,
    l2_error_2d,
    h1_semi_error_2d,
    quad_q1_rule,
    tri_p1_rule,
)
from beam_solution import BeamSolution2D

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


def _boundary_edges(nx, ny):
    """Return (bottom, top, left, right) edge lists for a structured nx×ny grid."""
    bottom = [(i, i + 1)                                    for i in range(nx - 1)]
    top    = [((ny - 1) * nx + i, (ny - 1) * nx + i + 1)    for i in range(nx - 1)]
    left   = [(j * nx, (j + 1) * nx)                        for j in range(ny - 1)]
    right  = [(j * nx + (nx - 1), (j + 1) * nx + (nx - 1))  for j in range(ny - 1)]
    return bottom, top, left, right


# ---------------------------------------------------------------------------
# Element strategies
#
# Each element type is a thin class declaring:
#   LABEL              : human-readable identifier (used in plot legends)
#   ELEMENT_RULE       : an element rule (e.g. quad_q1_rule(2))
#   add_topology(rootNode, Beam, nx, ny, L) : wire the SOFA topology and
#       return the topology container that holds the connectivity.
#       Uses RegularGridTopology under a sibling `Grid` child.
#   read_connectivity(topology) : extract the connectivity array from the
#       SOFA topology container, post-init.
#   to_triangles(conn) : triangulation for matplotlib tricontourf.
#
# Assembly / error-norm logic lives once on _ElementBase. Connectivity is
# supplied by the caller (the NodalForceAssembler controller reads it from
# SOFA after init); the element class never holds its own copy.
# ---------------------------------------------------------------------------

class _ElementBase:
    @classmethod
    def compute_nodal_forces(cls, nodes_2d, conn, mms, L, E, nu, nx, ny, dim):
        xy = nodes_2d[:, :2]

        F = assemble_nodal_forces_2d(
            lambda x, y: mms.source(x, y, E, nu, L, dim),
            xy, conn, cls._source_rule(mms))

        bottom, top, left, right = _boundary_edges(nx, ny)
        sides = [(bottom, 0.0, -1.0),
                 (top,    0.0, +1.0),
                 (left,  -1.0,  0.0),
                 (right, +1.0,  0.0)]
        for edges, nrm_x, nrm_y in sides:
            F += assemble_traction_2d(
                lambda x, y, nx=nrm_x, ny=nrm_y:
                    mms.traction(x, y, nx, ny, E, nu, L, dim),
                xy, edges)
        return F

    @classmethod
    def compute_l2(cls, sol, mms, L):
        return l2_error_2d(
            sol.nodes, sol.conn, np.column_stack([sol.ux, sol.uy]),
            lambda x, y: mms.u_ex(x, y, L),
            cls.ELEMENT_RULE)

    @classmethod
    def compute_h1(cls, sol, mms, L):
        return h1_semi_error_2d(
            sol.nodes, sol.conn, np.column_stack([sol.ux, sol.uy]),
            lambda x, y: mms.grad_u_ex(x, y, L),
            cls.ELEMENT_RULE)


class _QuadElement(_ElementBase):
    LABEL        = "Q1 quad"
    ELEMENT_RULE = staticmethod(quad_q1_rule(2))   # used for L²/H¹ error norms

    @staticmethod
    def _source_rule(mms):
        rule = mms.source_quadrature_quad
        if rule is None:
            raise ValueError(
                f"{type(mms).__name__}.source_quadrature_quad must be set")
        return rule

    @staticmethod
    def add_topology(Beam):
        topology = Beam.addObject("QuadSetTopologyContainer", name="topology",
                                  quads="@../Grid/grid.quads",
                                  position="@../Grid/grid.position")
        Beam.addObject("QuadSetTopologyModifier")
        return topology

    @staticmethod
    def read_connectivity(topology):
        return topology.quads.array().copy()

    @staticmethod
    def to_triangles(conn):
        tris = []
        for q in conn:
            tris.append([q[0], q[1], q[2]])
            tris.append([q[0], q[2], q[3]])
        return np.array(tris)


class _TriElement(_ElementBase):
    LABEL        = "P1 tri"
    ELEMENT_RULE = staticmethod(tri_p1_rule(3))    # used for L²/H¹ error norms

    @staticmethod
    def _source_rule(mms):
        rule = mms.source_quadrature_tri
        if rule is None:
            raise ValueError(
                f"{type(mms).__name__}.source_quadrature_tri must be set")
        return rule

    @staticmethod
    def add_topology(Beam):
        topology = Beam.addObject("TriangleSetTopologyContainer", name="topology")
        Beam.addObject("Quad2TriangleTopologicalMapping",
                       input="@../Grid/grid", output="@topology")
        Beam.addObject("TriangleSetTopologyModifier")
        return topology

    @staticmethod
    def read_connectivity(topology):
        return topology.triangles.array().copy()

    @staticmethod
    def to_triangles(conn):
        return conn


# ---------------------------------------------------------------------------
# Force assembly controller
#
# SOFA topology components (RegularGridTopology, QuadSetTopologyContainer,
# Quad2TriangleTopologicalMapping, ...) are only populated after init runs,
# not during Python scene-build time. This controller defers nodal-force
# assembly to onSimulationInitDoneEvent, where it reads positions off the
# MechanicalObject and connectivity off the SOFA topology, then fills the
# ConstantForceField that was added with placeholder zeros.
#
# Mirrors the 1D pattern in `bar.py:BodyForceAssembler`.
# ---------------------------------------------------------------------------

class NodalForceAssembler(Sofa.Core.Controller):
    def __init__(self, dofs, topology, force_field, element, mms,
                 L, E, nu, nx, ny, dim, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.dofs        = dofs
        self.topology    = topology
        self.force_field = force_field
        self.element     = element
        self.mms         = mms
        self.L, self.E, self.nu     = L, E, nu
        self.nx, self.ny, self.dim  = nx, ny, dim

    def onSimulationInitDoneEvent(self, event):
        nodes = self.dofs.rest_position.array().copy()
        conn  = self.element.read_connectivity(self.topology)
        F_xy  = self.element.compute_nodal_forces(
            nodes, conn, self.mms,
            self.L, self.E, self.nu, self.nx, self.ny, self.dim)
        if self.dim == "3d":
            F_full = np.hstack([F_xy, np.zeros((len(F_xy), 1))])
        else:
            F_full = F_xy
        with self.force_field.forces.writeableArray() as forces:
            forces[:] = F_full


# ---------------------------------------------------------------------------
# SOFA scene
# ---------------------------------------------------------------------------

def build_beam_scene(rootNode, mms, element, L=1.0, E=1e6, nu=0.3,
                     nx=10, ny=10, with_visual=True, dim="2d"):
    """Build a SOFA scene for `mms` on the 2D `element` strategy.

    dim : "2d" → Vec2d template / plane stress
          "3d" → Vec3d template / plane strain (z coordinate fixed at 0)
    Returns (dofs, topology). Nodes and connectivity become available
    after `Sofa.Simulation.init(root)` runs, via
    `dofs.rest_position.array()` and `element.read_connectivity(topology)`.
    """
    tmpl = _dim_template(dim)

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

    # Nodes still computed in Python so MechanicalObject's Vec2d/Vec3d template
    # receives the right per-point dimension; SOFA-side positions on the
    # topology container (via @../Grid/grid.position) are kept in sync by the
    # regular-grid construction.
    nodes_2d = get_nodes_2d(L, nx, ny, dim=dim)

    # Grid must be added *before* Beam: SOFA inits children in insertion
    # order, and Beam's topology container resolves `@../Grid/grid.position`
    # during its own init — so Grid has to be initialised first.
    Grid = rootNode.addChild("Grid")
    Grid.addObject("RegularGridTopology", name="grid",
                   nx=nx, ny=ny, nz=1,
                   min=[0.0, 0.0, 0.0], max=[L, L, 0.0])

    Beam = rootNode.addChild("Beam")
    Beam.addObject("StaticSolver", name="staticSolver", printLog=False)
    Beam.addObject("NewtonRaphsonSolver", name="newtonSolver",
                   maxNbIterationsNewton=1,
                   absoluteResidualStoppingThreshold=1e-10,
                   printLog=False)
    Beam.addObject("SparseLDLSolver", name="linearSolver",
                   template="CompressedRowSparseMatrixd")

    dofs = Beam.addObject("MechanicalObject", name="dofs", template=tmpl,
                          position=nodes_2d.tolist(),
                          showObject=with_visual, showObjectScale=0.005 * L)

    topology = element.add_topology(Beam)

    Beam.addObject("LinearSmallStrainFEMForceField", name="FEM", template=tmpl,
                   youngModulus=E, poissonRatio=nu, topology="@topology")

    mms.apply_bcs(Beam, nodes_2d, L, dim)

    # Placeholder force field filled in by the controller after init.
    n_nodes = len(nodes_2d)
    n_comp  = 3 if dim == "3d" else 2
    force_field = Beam.addObject("ConstantForceField", name="MMS_forces",
                                  template=tmpl,
                                  indices=list(range(n_nodes)),
                                  forces=[[0.0] * n_comp] * n_nodes)

    Beam.addObject(NodalForceAssembler(
        dofs=dofs, topology=topology, force_field=force_field,
        element=element, mms=mms,
        L=L, E=E, nu=nu, nx=nx, ny=ny, dim=dim,
        name="nodalForceAssembler"))

    return dofs, topology


# ─────────────────────────────────────────────────────────────────────────────
# Simulation runner
# ─────────────────────────────────────────────────────────────────────────────

def solve_beam(elem, mms, L, E, nu, nx, ny, dim="2d"):
    """Build, init, and run one static step. Returns a BeamSolution2D snapshot."""
    root = Sofa.Core.Node("root")
    dofs, topology = build_beam_scene(
        root, mms, elem, L=L, E=E, nu=nu, nx=nx, ny=ny, with_visual=False, dim=dim
    )
    Sofa.Simulation.init(root)
    # Read topology back from SOFA now that init has populated it.
    nodes_2d = dofs.rest_position.array().copy()
    conn     = elem.read_connectivity(topology)
    pos0     = dofs.position.array().copy()
    Sofa.Simulation.animate(root, root.dt.value)
    pos1     = dofs.position.array().copy()
    Sofa.Simulation.unload(root)
    # Extract x and y displacement columns (works for both Vec2d and Vec3d)
    ux = pos1[:, 0] - pos0[:, 0]
    uy = pos1[:, 1] - pos0[:, 1]
    return BeamSolution2D(nodes=nodes_2d, conn=conn, ux=ux, uy=uy)


# ---------------------------------------------------------------------------
# Output helpers (mirror 1D bar.py)
# ---------------------------------------------------------------------------

def write_solution_table(stem, x, y, ux_h, uy_h, u_ex, error_dict):
    """
    Write per-node solution table and error summary to results/<stem>.txt.

    u_ex       : callable (x, y) -> (ux_ex, uy_ex)
    error_dict : ordered dict of {label: value} for error summary lines
    """
    os.makedirs(RESULTS_DIR, exist_ok=True)
    path = os.path.join(RESULTS_DIR, f"{stem}.txt")
    with open(path, "w") as f:
        f.write(f"{'x':>10} | {'y':>10} | {'ux_h':>15} | {'uy_h':>15} | "
                f"{'ux_ex':>15} | {'uy_ex':>15} | {'err_x':>15} | {'err_y':>15}\n")
        f.write("-" * 124 + "\n")
        for xi, yi, uxi, uyi in zip(x, y, ux_h, uy_h):
            uxe, uye = u_ex(xi, yi)
            f.write(f"{xi:10.4f} | {yi:10.4f} | {uxi:15.6e} | {uyi:15.6e} | "
                    f"{uxe:15.6e} | {uye:15.6e} | "
                    f"{abs(uxi - uxe):15.6e} | {abs(uyi - uye):15.6e}\n")
        f.write("\n")
        for label, val in error_dict.items():
            f.write(f"{label:12s} = {val:.6e}\n")


def plot_solution_profile(stem, sol, mms, L, nx, ny, label, dim, hyp, nu, l2, h1):
    """Save 1-D mid-height profile (ux, uy vs x) to results/<stem>.png."""
    os.makedirs(RESULTS_DIR, exist_ok=True)
    xy    = sol.nodes[:, :2]
    mid_j = (ny - 1) // 2
    sl    = slice(mid_j * nx, mid_j * nx + nx)
    yc    = xy[mid_j * nx, 1]
    xf    = np.linspace(0, L, 300)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    for ax, u_sofa, u_fn, lbl, fmt in zip(
        axes,
        [sol.ux[sl], sol.uy[sl]],
        [lambda x: mms.u_ex(x, yc, L)[0], lambda x: mms.u_ex(x, yc, L)[1]],
        [r"$u_x$", r"$u_y$"],
        ["o-", "s-"],
    ):
        ax.plot(xy[sl, 0], u_sofa, fmt, color="tab:green",
                label=f"SOFA {label} [{dim}]", ms=5)
        ax.plot(xf, u_fn(xf), "--", color="tab:blue", label="MMS exact")
        ax.set_xlabel("x"); ax.set_ylabel(lbl)
        ax.legend(); ax.grid(True, alpha=0.3)
    fig.suptitle(f"{mms.name} — {label} [{dim} / {hyp}]  "
                 f"nu={nu}  nx={nx}  |L2={l2:.2e}  H1={h1:.2e}")
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, f"{stem}.png"), dpi=150)
    plt.close(fig)


def plot_solution_fields(stem, sol, elem, mms, L, nx, label, dim, hyp, nu):
    """Save 2-D field colour maps (ux, uy, error) to results/<stem>.png."""
    os.makedirs(RESULTS_DIR, exist_ok=True)
    xy             = sol.nodes[:, :2]
    ux_ref, uy_ref = mms.u_ex(xy[:, 0], xy[:, 1], L)
    tris_plot      = elem.to_triangles(sol.conn)
    x, y           = xy[:, 0], xy[:, 1]

    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    for ax, data, title, cmap in [
        (axes[0,0], sol.ux,                r"$u_x$ SOFA",            "gray"),
        (axes[0,1], ux_ref,                r"$u_x$ MMS",             "gray"),
        (axes[0,2], np.abs(sol.ux-ux_ref), r"$|u_x - u_x^{MMS}|$",  "gray_r"),
        (axes[1,0], sol.uy,                r"$u_y$ SOFA",            "gray"),
        (axes[1,1], uy_ref,                r"$u_y$ MMS",             "gray"),
        (axes[1,2], np.abs(sol.uy-uy_ref), r"$|u_y - u_y^{MMS}|$",  "gray_r"),
    ]:
        tc = ax.tricontourf(x, y, tris_plot.tolist(), data, levels=20, cmap=cmap)
        ax.triplot(x, y, tris_plot.tolist(), "k-", lw=0.3, alpha=0.4)
        plt.colorbar(tc, ax=ax, shrink=0.8)
        ax.set_title(title); ax.set_aspect("equal")
        ax.set_xlabel("x"); ax.set_ylabel("y")
    fig.suptitle(f"Fields 2D — {label} [{dim} / {hyp}]  "
                 f"{mms.name}  nu={nu}  nx={nx}")
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, f"{stem}.png"), dpi=150)
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# Single-case driver
# ─────────────────────────────────────────────────────────────────────────────

def run_reference_scene(elem, mms):
    """Solve one MMS case at the reference mesh, write the solution table and plot.

    All parameters come from params.json (top-level + `reference` block).
    """
    cfg = load_params()
    ref = cfg["reference"]
    L, E    = cfg["length"], cfg["youngModulus"]
    nu, dim = ref["nu"], ref["dim"]
    nx = ny = ref["nx"]
    hyp     = "plane strain" if dim == "3d" else "plane stress"

    sol = solve_beam(elem, mms, L, E, nu, nx, ny, dim=dim)
    l2  = elem.compute_l2(sol, mms, L)
    h1  = elem.compute_h1(sol, mms, L)

    label = elem.LABEL
    tag   = label.replace(" ", "_")
    stem  = f"{mms.name}_{tag}_{dim}_nu{nu}_nx{nx}"

    xy = sol.nodes[:, :2]
    write_solution_table(f"solution_{stem}", xy[:, 0], xy[:, 1],
                         sol.ux, sol.uy,
                         lambda xi, yi: mms.u_ex(xi, yi, L),
                         {"L2": l2, "H1_semi": h1})
    plot_solution_profile(f"solution_{stem}", sol, mms, L, nx, ny,
                          label, dim, hyp, nu, l2, h1)
    plot_solution_fields (f"fields2D_{stem}", sol, elem, mms, L, nx,
                          label, dim, hyp, nu)


def case_scene(mms, element):
    """Return a `createScene(rootNode)` bound to this MMS and element type.

    All parameters come from params.json (top-level + `reference` block). Each
    case file exposes:
        createScene = case_scene(mms, element_quad)
    so that `runSofa cubic.py` loads the default scene.
    """
    def createScene(rootNode):
        cfg = load_params()
        ref = cfg["reference"]
        build_beam_scene(rootNode, mms, element,
                         L=cfg["length"], E=cfg["youngModulus"],
                         nu=ref["nu"], nx=ref["nx"], ny=ref["nx"],
                         with_visual=True, dim=ref["dim"])
        return rootNode
    return createScene


# ─────────────────────────────────────────────────────────────────────────────
# Element instances
# ─────────────────────────────────────────────────────────────────────────────

element_quad = _QuadElement()
element_tri  = _TriElement()

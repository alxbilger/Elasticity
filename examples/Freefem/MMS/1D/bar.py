"""1D bar problem: scene build, solve, post-processing for MMS validation."""

import json
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import Sofa
import Sofa.Core
import Sofa.Simulation

# Make the parent MMS/ directory importable so we can pull in fem.py.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from fem import (
    line_quadrature,
    assemble_nodal_forces_1d,
    l2_error_1d,
    h1_semi_error_1d,
    L2_QUADRATURE_1D,
    H1_QUADRATURE_1D,
)
from bar_solution import BarSolution1D
from output import write_solution_table
from scene import NodalForceAssembler

RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")


# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

def load_params(path=None):
    if path is None:
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "params.json")
    with open(path) as f:
        cfg = json.load(f)
    # Non-dimensional effective Young's modulus on x ∈ [0,1]: E_eff = E / L.
    cfg["E_eff"] = cfg["youngModulus"] / cfg["length"]
    return cfg


# ---------------------------------------------------------------------------
# SOFA runner
# ---------------------------------------------------------------------------

def _bar_force_compute(f_body, quadrature):
    """Return a `(nodes, topology) -> (nx, 1) force array` for the assembler."""
    def compute(nodes, topology):
        flat  = nodes.copy().flatten()
        edges = topology.edges.array().copy()
        return assemble_nodal_forces_1d(f_body, flat, edges, quadrature).reshape(-1, 1)
    return compute


def build_bar_scene(root, mms, E_eff, nx):
    """Populate root with a static 1D bar scene on the non-dimensional domain [0,1].

    BodyForce is assembled after init by BodyForceAssembler. BCs: Dirichlet at
    x=0, Neumann `mms.traction_bc(E_eff)` at x=1.
    """
    root.addObject('RequiredPlugin', pluginName=[
        "Elasticity",
        "Sofa.Component.Constraint.Projective",
        "Sofa.Component.LinearSolver.Direct",
        "Sofa.Component.MechanicalLoad",
        "Sofa.Component.ODESolver.Backward",
        "Sofa.Component.StateContainer",
        "Sofa.Component.Topology.Container.Grid",
        "Sofa.Component.Topology.Container.Dynamic",
        "Sofa.Component.Visual",
    ])
    root.addObject('DefaultAnimationLoop')
    root.addObject('VisualStyle',
                   displayFlags="showBehaviorModels showForceFields")

    Grid = root.addChild('Grid')
    Grid.addObject('RegularGridTopology',
                  name="grid",
                  nx=nx, ny=1, nz=1,
                  min=[0., 0., 0.],
                  max=[1., 0., 0.])

    Bar = root.addChild('Bar')
    Bar.addObject('NewtonRaphsonSolver',
                  name="newtonSolver",
                  maxNbIterationsNewton=10,
                  absoluteResidualStoppingThreshold=1e-10,
                  printLog=False)
    Bar.addObject('SparseLDLSolver',
                  name="linearSolver",
                  template="CompressedRowSparseMatrixd")
    Bar.addObject('StaticSolver',
                  name="staticSolver",
                  newtonSolver="@newtonSolver",
                  linearSolver="@linearSolver")

    Bar.addObject('EdgeSetTopologyContainer', name="topology",
                  edges="@../Grid/grid.edges",
                  position="@../Grid/grid.position")

    dofs = Bar.addObject('MechanicalObject',
                         name="dofs",
                         template="Vec1d")

    Bar.addObject('LinearSmallStrainFEMForceField',
                  name="FEM",
                  template="Vec1d",
                  youngModulus=E_eff,
                  topology="@topology")

    body_force = Bar.addObject('ConstantForceField',
                               name="BodyForce",
                               indices=list(range(nx)),
                               forces=[[0.0]] * nx)

    Bar.addObject(NodalForceAssembler(
        dofs=dofs,
        topology=Bar.topology,
        force_field=body_force,
        compute_forces=_bar_force_compute(
            f_body=lambda xi: mms.source(xi, E_eff),
            quadrature=mms.source_quadrature),
        name="bodyForceAssembler"))

    mms.apply_bcs(Bar, E_eff, nx)


def case_scene(mms):
    """Return a `createScene(rootNode)` bound to this MMS case."""
    def createScene(rootNode):
        cfg = load_params()
        build_bar_scene(rootNode, mms, cfg["E_eff"], cfg["reference"]["nx"])
        return rootNode
    return createScene


def solve_bar(mms, E_eff, nx):
    """Build, run one static step, and return a BarSolution1D snapshot."""
    root = Sofa.Core.Node("root")
    build_bar_scene(root, mms, E_eff, nx)
    Sofa.Simulation.init(root)
    Sofa.Simulation.animate(root, root.dt.value)
    Bar   = root.Bar
    x0    = Bar.dofs.rest_position.array().copy().flatten()
    edges = Bar.topology.edges.array().copy()
    u_h   = Bar.dofs.position.array().copy().flatten() - x0
    Sofa.Simulation.unload(root)
    return BarSolution1D(x0=x0, edges=edges, u_h=u_h)


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def plot_solution(case, x, u_h, u_ex, label_ex):
    """Save solution comparison plot to results/<case>_solution.png.

    u_ex : callable x -> float
    """
    os.makedirs(RESULTS_DIR, exist_ok=True)
    x_fine = np.linspace(x[0], x[-1], 300)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(x_fine, u_ex(x_fine), "k--", label=label_ex, linewidth=2)
    ax.plot(x,      u_h,        "bo",  label="SOFA",   markersize=6)
    ax.set_xlabel("x")
    ax.set_ylabel("u(x)")
    ax.set_title(f"MMS {case} - SOFA vs exact")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, f"{case}_solution.png"), dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Single-case driver
# ---------------------------------------------------------------------------

def run_reference_scene(mms):
    """Solve one MMS case at the reference mesh, write the solution table and plot."""
    cfg = load_params()
    sol = solve_bar(mms, cfg["E_eff"], cfg["reference"]["nx"])
    l2  = l2_error_1d(sol.x0, sol.edges, sol.u_h, mms.u_ex, L2_QUADRATURE_1D)
    h1  = h1_semi_error_1d(sol.x0, sol.edges, sol.u_h, mms.du_ex, H1_QUADRATURE_1D)
    write_solution_table(f"solution_{mms.name}", sol.x0, sol.u_h, mms.u_ex,
                         RESULTS_DIR, {"L2": l2, "H1_semi": h1})
    plot_solution(mms.name, sol.x0, sol.u_h, mms.u_ex, mms.plot_label)

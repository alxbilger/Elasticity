"""Shared utilities for 1D MMS validation."""

import json
import os
import numpy as np
import matplotlib.pyplot as plt
import Sofa
import Sofa.Core
import Sofa.Simulation

RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")


# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

def load_params(path=None):
    if path is None:
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "params.json")
    with open(path) as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Quadrature rules
# Approximation: integral_{x_i}^{x_{i+1}} g(x) dx
#                ~= (h/2) * sum_k( w_k * g(x_k) )
# where x_k = (x_i + x_{i+1})/2 + (h/2)*xi_k
# ---------------------------------------------------------------------------

def midpoint_quadrature(g, x1, x2):
    """1-point Gauss rule: n_g=1, xi=[0], w=[2]."""
    h   = x2 - x1
    x_k = 0.5 * (x1 + x2)
    return (h / 2.0) * 2.0 * g(x_k)


def gauss2_quadrature(g, x1, x2):
    """2-point Gauss rule: n_g=2, xi=[-1/sqrt(3), +1/sqrt(3)], w=[1, 1]."""
    h     = x2 - x1
    x_mid = 0.5 * (x1 + x2)
    xi    = np.array([-1.0 / np.sqrt(3.0), 1.0 / np.sqrt(3.0)])
    x_k   = x_mid + (h / 2.0) * xi
    return (h / 2.0) * (g(x_k[0]) + g(x_k[1]))


# ---------------------------------------------------------------------------
# FEM assembly
# ---------------------------------------------------------------------------

def assemble_nodal_forces(f_body, nodes, quadrature):
    """
    Assemble the consistent nodal force vector F_i = integral f_body(x) phi_i(x) dx.

    f_body    : callable x -> float
    nodes     : 1-D array of node coordinates
    quadrature: midpoint_quadrature or gauss2_quadrature
    """
    forces = np.zeros(len(nodes))
    for i in range(len(nodes) - 1):
        x1, x2 = nodes[i], nodes[i + 1]
        h = x2 - x1
        forces[i]     += quadrature(lambda x, x1=x1, x2=x2, h=h: f_body(x) * (x2 - x) / h, x1, x2)
        forces[i + 1] += quadrature(lambda x, x1=x1, x2=x2, h=h: f_body(x) * (x - x1) / h, x1, x2)
    return forces


# ---------------------------------------------------------------------------
# Error norms
# ---------------------------------------------------------------------------

def l2_error(nodes, u_h, u_ex, quadrature):
    """L2 error norm: sqrt( integral (u_h - u_ex)^2 dx ) over the mesh."""
    total = 0.0
    for i in range(len(nodes) - 1):
        x1, x2 = nodes[i], nodes[i + 1]
        h = x2 - x1
        u_a, u_b = u_h[i], u_h[i + 1]
        u_interp = lambda x, x1=x1, h=h, u_a=u_a, u_b=u_b: u_a + (u_b - u_a) * (x - x1) / h
        total += quadrature(lambda x: (u_interp(x) - u_ex(x)) ** 2, x1, x2)
    return np.sqrt(total)


def h1_semi_error(nodes, u_h, du_ex, quadrature):
    """H1 semi-norm error: sqrt( integral (du_h - du_ex)^2 dx ) over the mesh."""
    total = 0.0
    for i in range(len(nodes) - 1):
        x1, x2 = nodes[i], nodes[i + 1]
        h = x2 - x1
        # du_h = sum_a u_a * dphi_a/dx, with dphi_a/dx = -1/h, dphi_b/dx = +1/h
        du_h = u_h[i] * (-1.0 / h) + u_h[i + 1] * (1.0 / h)
        total += quadrature(lambda x, du_h=du_h: (du_h - du_ex(x)) ** 2, x1, x2)
    return np.sqrt(total)


# ---------------------------------------------------------------------------
# SOFA runner
# ---------------------------------------------------------------------------


def build_bar_scene(root, young_modulus, nx, nodal_forces, apply_bcs, L=1.0):
    """
    Populate root with a static 1D bar scene on the non-dimensional domain [0, L].
    Uses RegularGridTopology for mesh definition.

    Parameters
    ----------
    root         : SOFA root node
    young_modulus: float
    nx           : int, number of nodes along the bar
    nodal_forces : array of length nx, assembled body force vector
    apply_bcs    : callable(Bar, nx) that adds all boundary conditions
                   (FixedProjectiveConstraint + Neumann ConstantForceField)
    L            : float, bar length (default 1.0, non-dimensional)

    Returns the dofs MechanicalObject.
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

    root.gravity.value = [0, 0, 0]
    root.dt.value = 1.0

    Bar = root.addChild('Bar')

    #  Regular grid Topology 
    Bar.addObject('RegularGridTopology',
                  name="grid",
                  nx=nx, ny=1, nz=1,
                  min=[0., 0., 0.],
                  max=[L, 0., 0.])

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

    dofs = Bar.addObject('MechanicalObject',
                         name="dofs",
                         template="Vec1d")

    Bar.addObject('EdgeSetTopologyContainer',
                  name="topology",
                  src="@grid")

    Bar.addObject('LinearSmallStrainFEMForceField',
                  name="FEM",
                  template="Vec1d",
                  youngModulus=young_modulus,
                  topology="@topology")

    Bar.addObject('ConstantForceField',
                  name="BodyForce",
                  indices=list(range(nx)),
                  forces=nodal_forces)

    apply_bcs(Bar, nx)

    return dofs


def run_bar_mms(young_modulus, nx, nodal_forces, apply_bcs):
    """Build, run one static step, and return (node_positions, displacements)."""
    root = Sofa.Core.Node("root")
    dofs = build_bar_scene(root, young_modulus, nx, nodal_forces, apply_bcs)
    Sofa.Simulation.init(root)
    x0 = dofs.position.array().copy().flatten()
    Sofa.Simulation.animate(root, root.dt.value)
    u  = dofs.position.array().flatten() - x0
    Sofa.Simulation.unload(root)
    return x0, u


# ---------------------------------------------------------------------------
# Convergence study
# ---------------------------------------------------------------------------

def convergence_study(case, nx_list, run_fn, error_fns, ref_slopes):
    """
    Run a mesh refinement convergence study and write results.

    run_fn    : callable(nx) -> (x0, u_h)
    error_fns : dict { label: callable(x0, u_h) -> float }
    ref_slopes: dict { label: (error_key, slope) }
                where error_key is a key in error_fns, used to anchor the
                O(h^p) reference line in the convergence plot
    """
    hs     = []
    errors = {label: [] for label in error_fns}
    rows   = []

    for k, nx_k in enumerate(nx_list):
        h_k     = 1.0 / (nx_k - 1)
        x0, u_h = run_fn(nx_k)
        row     = {"nx": nx_k, "h": h_k}

        for label, err_fn in error_fns.items():
            e_k  = err_fn(x0, u_h)
            rate = (f"{np.log(e_k / errors[label][-1]) / np.log(h_k / hs[-1]):.2f}"
                    if k > 0 else "")
            errors[label].append(e_k)
            row[label]            = e_k
            row[f"rate_{label}"]  = rate

        hs.append(h_k)
        rows.append(row)

    write_convergence_table(case, rows)
    plot_convergence(case, hs, errors,
                     {label: (errors[err_key], slope)
                      for label, (err_key, slope) in ref_slopes.items()})


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def write_solution_table(case, x, u_h, u_ex, error_dict):
    """
    Write per-node solution table and error summary to results/<case>_solution.txt.

    u_ex : callable x -> float
    error_dict : ordered dict of {label: value} for error summary lines
    """
    os.makedirs(RESULTS_DIR, exist_ok=True)
    path = os.path.join(RESULTS_DIR, f"{case}_solution.txt")
    with open(path, "w") as f:
        f.write(f"{'x':>10} | {'u_h':>15} | {'u_ex':>15} | {'error':>15}\n")
        f.write("-" * 60 + "\n")
        for xi, ui in zip(x, u_h):
            ue  = u_ex(xi)
            f.write(f"{xi:10.4f} | {ui:15.6e} | {ue:15.6e} | {abs(ui - ue):15.6e}\n")
        f.write("\n")
        for label, val in error_dict.items():
            f.write(f"{label:12s} = {val:.6e}\n")


def write_convergence_table(case, rows):
    """
    Write convergence table to results/<case>_convergence.txt.

    rows : list of dicts with keys 'nx', 'h', and one key per error column.
           Rate columns are strings (empty for the first row).
    """
    os.makedirs(RESULTS_DIR, exist_ok=True)
    path = os.path.join(RESULTS_DIR, f"{case}_convergence.txt")
    err_keys = [k for k in rows[0] if k not in ("nx", "h")]
    header   = f"{'nx':>6} | {'h':>10}" + "".join(f" | {k:>16}" for k in err_keys)
    with open(path, "w") as f:
        f.write(header + "\n")
        f.write("-" * len(header) + "\n")
        for row in rows:
            line = f"{row['nx']:6d} | {row['h']:10.4f}"
            for k in err_keys:
                v = row[k]
                line += f" | {v:16.6e}" if isinstance(v, float) else f" | {v:>16}"
            f.write(line + "\n")


def plot_solution(case, x, u_h, u_ex, label_ex):
    """Save solution comparison plot to results/<case>_solution.png.

    u_ex : callable x -> float
    """
    os.makedirs(RESULTS_DIR, exist_ok=True)
    x_fine = np.linspace(x[0], x[-1], 300)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(x,      u_h,        "bo-", label="SOFA",   markersize=6, linewidth=2)
    ax.plot(x_fine, u_ex(x_fine), "r--", label=label_ex, linewidth=2)
    ax.set_xlabel("x")
    ax.set_ylabel("u(x)")
    ax.set_title(f"MMS {case}: SOFA vs exact")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, f"{case}_solution.png"), dpi=150)
    plt.close(fig)


def plot_convergence(case, hs, error_series, reference_slopes):
    """
    Save log-log convergence plot to results/<case>_convergence.png.

    error_series    : dict { label: array-like of errors }
    reference_slopes: dict { label: (errors_ref, slope) }
    """
    os.makedirs(RESULTS_DIR, exist_ok=True)
    h_arr   = np.array(hs)
    h_ref   = np.array([h_arr[0], h_arr[-1]])
    markers = ["bo-", "rs--", "g^:"]
    fig, ax = plt.subplots(figsize=(8, 5))
    for (label, errors), marker in zip(error_series.items(), markers):
        ax.loglog(h_arr, errors, marker, label=label, linewidth=2, markersize=7)
    for label, (ref_errors, slope) in reference_slopes.items():
        e_last = np.array(ref_errors)[-1]
        ax.loglog(h_ref, e_last * (h_ref / h_arr[-1]) ** slope,
                  "k--", linewidth=1.5, label=label)
    ax.set_xlabel("h")
    ax.set_ylabel("Error")
    ax.set_title(f"Convergence — {case}")
    ax.legend()
    ax.grid(True, alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, f"{case}_convergence.png"), dpi=150)
    plt.close(fig)

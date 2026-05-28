"""Run the mesh-refinement convergence study for every 2D MMS case.

Loops cases × elements × dim × nu × nx and writes per-(case, element, dim, nu)
text tables and convergence plots into the shared `results/` directory.
"""

from cubic          import mms as cubic_mms
from trigonometric  import mms as trig_mms
from incompressible import mms as incomp_mms

from beam import (
    RESULTS_DIR,
    load_params,
    solve_beam,
    element_quad,
    element_tri,
)
from convergence import run_convergence_series
from output      import plot_convergence


def convergence_study(elem_specs, mms, L, E, nu, nx_values, dim="2d"):
    """
    Run a convergence series for each element type in elem_specs, write a
    per-(element) text table, and one shared plot with L²/H¹ for every
    element on the same axes.

    elem_specs : list of dicts with keys 'elem', 'label', 'l2_style', 'h1_style'
    dim        : "2d" (plane stress) or "3d" (plane strain)
    """
    hyp = "plane strain" if dim == "3d" else "plane stress"
    print(f"\n  PoissonRatio = {nu}  ({hyp})", flush=True)

    plot_series, hs_ref = [], None
    for spec in elem_specs:
        elem, label = spec["elem"], spec["label"]
        tag         = label.replace(" ", "_")
        stem        = f"convergence_{mms.name}_{tag}_{dim}_nu{nu}"

        hs, errors = run_convergence_series(
            nx_values  = nx_values,
            run_fn     = lambda nx, _e=elem: solve_beam(
                _e, mms, L, E, nu, nx, nx, dim=dim),
            h_fn       = lambda nx: L / (nx - 1),
            error_fns  = {
                "L2": lambda sol, _e=elem: _e.compute_l2(sol, mms, L),
                "H1": lambda sol, _e=elem: _e.compute_h1(sol, mms, L),
            },
            banner     = f"-- {label}  [{dim} / {hyp}]  {mms.name}  nu={nu} --",
            results_dir = RESULTS_DIR,
            table_stem  = stem,
        )

        plot_series.append({"label": f"{label} L²",
                            "errors": errors["L2"], "style": spec["l2_style"]})
        plot_series.append({"label": f"{label} H¹",
                            "errors": errors["H1"], "style": spec["h1_style"]})
        hs_ref = hs

    title = f"Convergence — {mms.name} [{dim} / {hyp}]  nu={nu}"
    plot_convergence(f"convergence_{mms.name}_{dim}_nu{nu}",
                     hs_ref, plot_series, title=title, results_dir=RESULTS_DIR)


if __name__ == "__main__":
    cfg  = load_params()
    L    = cfg["length"]
    E    = cfg["youngModulus"]
    conv = cfg["convergence"]

    specs = [
        {"elem": element_quad, "label": "Q1 quad",
         "l2_style": "bo-",  "h1_style": "rs--"},
        {"elem": element_tri,  "label": "P1 tri",
         "l2_style": "b^-",  "h1_style": "rD--"},
    ]

    # dim="2d" → Vec2d / plane stress;  dim="3d" → Vec3d / plane strain
    for mms in (cubic_mms, trig_mms, incomp_mms):
        nx_vals = conv["nx_values"][mms.name]
        print(f"\n== {mms.name} ==")
        for DIM in conv["dim_values"]:
            for nu in conv["nu_values"]:
                convergence_study(specs, mms, L, E, nu, nx_vals, dim=DIM)

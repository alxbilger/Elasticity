"""Run the mesh-refinement convergence study for every 3D MMS case.

Loops cases × nu × nx and writes per-(case, nu) text tables and convergence
plots into the shared `results/` directory. Mirrors the 2D driver minus the
plane-stress / plane-strain `dim` axis (3D has a single constitutive branch).
"""

from sinus_neumann import mms as sinus_neumann_mms

from solid import (
    RESULTS_DIR,
    load_params,
    solve_solid,
    element_hex,
)
from convergence import run_convergence_series
from output      import plot_convergence


def convergence_study(elem_specs, mms, L, E, nu, nx_values):
    """
    Run a convergence series for each element type in elem_specs, write a
    per-(element) text table, and one shared plot with L²/H¹ for every
    element on the same axes.

    elem_specs : list of dicts with keys 'elem', 'label', 'l2_style', 'h1_style'
    """
    print(f"\n  PoissonRatio = {nu}", flush=True)

    plot_series, hs_ref = [], None
    for spec in elem_specs:
        elem, label = spec["elem"], spec["label"]
        tag         = label.replace(" ", "_")
        stem        = f"convergence_{mms.name}_{tag}_nu{nu}"

        hs, errors = run_convergence_series(
            nx_values  = nx_values,
            run_fn     = lambda nx, _e=elem: solve_solid(
                _e, mms, L, E, nu, nx, nx, nx),
            h_fn       = lambda nx: L / (nx - 1),
            error_fns  = {
                "L2": lambda sol, _e=elem: _e.compute_l2(sol, mms, L),
                "H1": lambda sol, _e=elem: _e.compute_h1(sol, mms, L),
            },
            banner     = f"-- {label}  {mms.name}  nu={nu} --",
            results_dir = RESULTS_DIR,
            table_stem  = stem,
        )

        plot_series.append({"label": f"{label} L²",
                            "errors": errors["L2"], "style": spec["l2_style"]})
        plot_series.append({"label": f"{label} H¹",
                            "errors": errors["H1"], "style": spec["h1_style"]})
        hs_ref = hs

    title = f"Convergence — {mms.name}  nu={nu}"
    plot_convergence(f"convergence_{mms.name}_nu{nu}",
                     hs_ref, plot_series, title=title, results_dir=RESULTS_DIR)


if __name__ == "__main__":
    cfg  = load_params()
    L    = cfg["length"]
    E    = cfg["youngModulus"]
    conv = cfg["convergence"]

    specs = [
        {"elem": element_hex, "label": "Q1 hex",
         "l2_style": "bo-", "h1_style": "rs--"},
    ]

    for mms in (sinus_neumann_mms,):
        nx_vals = conv["nx_values"][mms.name]
        print(f"\n== {mms.name} ==")
        for nu in conv["nu_values"]:
            convergence_study(specs, mms, L, E, nu, nx_vals)

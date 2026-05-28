"""Run the mesh-refinement convergence study for every 1D MMS case."""

from quadratic   import mms as quadratic_mms
from cubic       import mms as cubic_mms
from sinusoidal  import mms as sinusoidal_mms
from exponential import mms as exponential_mms

from bar import (
    RESULTS_DIR,
    load_params,
    solve_bar,
    l2_error_1d,
    h1_semi_error_1d,
    L2_QUADRATURE_1D,
    H1_QUADRATURE_1D,
)
from convergence import run_convergence_series
from output      import plot_convergence


if __name__ == "__main__":
    cfg = load_params()
    for mms in (quadratic_mms, cubic_mms, sinusoidal_mms, exponential_mms):
        stem = f"{mms.name}_convergence"
        hs, errors = run_convergence_series(
            nx_values  = cfg["convergence"]["nx_values"][mms.name],
            run_fn     = lambda nx, _m=mms: solve_bar(_m, cfg["E_eff"], nx),
            h_fn       = lambda nx: 1.0 / (nx - 1),
            error_fns  = {
                "L2": lambda sol, _m=mms: l2_error_1d(
                    sol.x0, sol.edges, sol.u_h, _m.u_ex, L2_QUADRATURE_1D),
                "H1": lambda sol, _m=mms: h1_semi_error_1d(
                    sol.x0, sol.edges, sol.u_h, _m.du_ex, H1_QUADRATURE_1D),
            },
            banner     = f"-- Convergence  {mms.name} --",
            results_dir = RESULTS_DIR,
            table_stem  = stem,
        )
        plot_series = [{"label": lbl, "errors": errors[lbl]} for lbl in errors]
        plot_convergence(stem, hs, plot_series,
                         title=f"Error Convergence — {mms.name} function",
                         results_dir=RESULTS_DIR)

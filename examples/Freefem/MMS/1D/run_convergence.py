"""Run the mesh-refinement convergence study for every 1D MMS case."""

import os
import numpy as np
import matplotlib.pyplot as plt

from quadratic  import mms as quadratic_mms
from cubic      import mms as cubic_mms
from sinusoidal import mms as sinusoidal_mms
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
from plot_utils import annotate_convergence_rates


def convergence_study(case, nx_list, run_fn, error_fns):
    """
    Run a mesh refinement convergence study and write results.

    run_fn    : callable(nx) -> BarSolution1D
    error_fns : dict { label: callable(x0, edges, u_h) -> float }
    """
    hs     = []
    errors = {label: [] for label in error_fns}
    rows   = []

    err_labels = list(error_fns.keys())
    hdr = f"{'nx':>5} | {'h':>10}" + "".join(
        f" | {lbl:>14} | {('rate_' + lbl):>7}" for lbl in err_labels)
    print(f"\n── Convergence  {case} ──\n{hdr}", flush=True)

    for k, nx_k in enumerate(nx_list):
        h_k = 1.0 / (nx_k - 1)
        sol = run_fn(nx_k)
        row = {"nx": nx_k, "h": h_k}

        for label, err_fn in error_fns.items():
            e_k  = err_fn(sol.x0, sol.edges, sol.u_h)
            rate = (f"{np.log(e_k / errors[label][-1]) / np.log(h_k / hs[-1]):.2f}"
                    if k > 0 else "")
            errors[label].append(e_k)
            row[label]            = e_k
            row[f"rate_{label}"]  = rate

        line = f"{nx_k:5d} | {h_k:10.4f}"
        for lbl in err_labels:
            line += f" | {row[lbl]:14.6e} | {row[f'rate_{lbl}']:>7}"
        print(line, flush=True)

        hs.append(h_k)
        rows.append(row)

    stem = f"{case}_convergence"
    write_convergence_table(stem, rows)
    plot_series = [{"label": lbl, "errors": errors[lbl]} for lbl in err_labels]
    plot_convergence(stem, hs, plot_series,
                     title=f"Error Convergence — {case} function")


def write_convergence_table(stem, rows):
    """
    Write convergence table to results/<stem>.txt.

    rows : list of dicts with keys 'nx', 'h', and one key per error column.
           Rate columns are strings (empty for the first row).
    """
    os.makedirs(RESULTS_DIR, exist_ok=True)
    path = os.path.join(RESULTS_DIR, f"{stem}.txt")
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


def plot_convergence(stem, hs, series, title, ylabel="Error"):
    """
    Save log-log convergence plot to results/<stem>.png.

    series : list of {"label", "errors", "style"?} dicts
    Per-segment convergence rates are annotated above each line segment.
    """
    os.makedirs(RESULTS_DIR, exist_ok=True)
    h_arr   = np.array(hs)
    default = ["bo-", "rs--", "g^:", "m^-"]
    fig, ax = plt.subplots(figsize=(8, 5))
    for i, s in enumerate(series):
        style = s.get("style", default[i % len(default)])
        e_arr = np.array(s["errors"])
        ax.loglog(h_arr, e_arr, style, label=s["label"], linewidth=2, markersize=7)
        annotate_convergence_rates(ax, h_arr, e_arr)
    ax.set_xlabel("h")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, f"{stem}.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    cfg = load_params()
    for mms in (quadratic_mms, cubic_mms, sinusoidal_mms, exponential_mms):
        convergence_study(mms.name, cfg["convergence"]["nx_values"][mms.name],
            run_fn    = lambda nx: solve_bar(mms, cfg["E_eff"], nx),
            error_fns = {"L2":          lambda x, e, u: l2_error_1d(x, e, u, mms.u_ex, L2_QUADRATURE_1D),
                         "H1 seminorm": lambda x, e, u: h1_semi_error_1d(x, e, u, mms.du_ex, H1_QUADRATURE_1D)})

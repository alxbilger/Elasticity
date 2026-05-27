"""Run the mesh-refinement convergence study for every 2D MMS case.

Loops cases × elements × dim × nu × nx and writes per-(case, element, dim, nu)
text tables and convergence plots into the shared `results/` directory.
"""

import os
import numpy as np
import matplotlib.pyplot as plt

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
from plot_utils import annotate_convergence_rates


def convergence_study(elem_specs, mms, L, E, nu, nx_values, dim="2d"):
    """
    Run convergence study for each element type in elem_specs, write a
    per-(element) text table, and one shared plot with L²/H¹ for every
    element on the same axes.

    elem_specs : list of dicts with keys 'elem', 'label', 'l2_style', 'h1_style'
    dim        : "2d" (plane stress) or "3d" (plane strain)
    """
    hyp = "plane strain" if dim == "3d" else "plane stress"
    print(f"\n  PoissonRatio = {nu}  ({hyp})", flush=True)
    hdr = (f"{'nx':>5} | {'h':>10} | {'L2':>14} | {'rate_L2':>7} "
           f"| {'H1':>14} | {'rate_H1':>7}")

    plot_series, hs_ref = [], None
    for spec in elem_specs:
        elem, label = spec["elem"], spec["label"]
        tag         = label.replace(" ", "_")
        stem        = f"convergence_{mms.name}_{tag}_{dim}_nu{nu}"

        print(f"\n── {label}  [{dim} / {hyp}]  {mms.name}  nu={nu} ──\n{hdr}",
              flush=True)

        rows, hs, l2s, h1s = [], [], [], []
        for k, nx in enumerate(nx_values):
            ny  = nx
            h   = L / (nx - 1)
            sol = solve_beam(elem, mms, L, E, nu, nx, ny, dim=dim)
            l2  = elem.compute_l2(sol, mms, L)
            h1  = elem.compute_h1(sol, mms, L)

            rate_l2 = (f"{np.log(l2 / l2s[-1]) / np.log(h / hs[-1]):.2f}"
                       if k > 0 else "")
            rate_h1 = (f"{np.log(h1 / h1s[-1]) / np.log(h / hs[-1]):.2f}"
                       if k > 0 else "")
            print(f"{nx:5d} | {h:10.4f} | {l2:14.6e} | {rate_l2:>7} "
                  f"| {h1:14.6e} | {rate_h1:>7}", flush=True)
            rows.append({"nx": nx, "h": h,
                         "L2": l2, "rate_L2": rate_l2,
                         "H1": h1, "rate_H1": rate_h1})
            hs.append(h); l2s.append(l2); h1s.append(h1)

        write_convergence_table(stem, rows)
        plot_series.append({"label": f"{label} L²",
                            "errors": l2s, "style": spec["l2_style"]})
        plot_series.append({"label": f"{label} H¹",
                            "errors": h1s, "style": spec["h1_style"]})
        hs_ref = hs

    title = f"Convergence — {mms.name} [{dim} / {hyp}]  nu={nu}"
    plot_convergence(f"convergence_{mms.name}_{dim}_nu{nu}",
                     hs_ref, plot_series, title=title)


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
        print(f"\n══ {mms.name} ══")
        for DIM in conv["dim_values"]:
            for nu in conv["nu_values"]:
                convergence_study(specs, mms, L, E, nu, nx_vals, dim=DIM)

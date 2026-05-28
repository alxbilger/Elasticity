"""Shared output helpers across MMS 1D/2D/3D drivers (tables + plots)."""

import os
import numpy as np
import matplotlib.pyplot as plt

from plot_utils import annotate_convergence_rates


_AXES = ("x", "y", "z")


def write_solution_table(stem, coords, u_h, u_ex, results_dir, error_dict):
    """
    Write per-node solution table and error summary to <results_dir>/<stem>.txt.

    coords     : (N,) or (N, dim) array of node coordinates
    u_h        : (N,) or (N, dim) array of nodal displacements
    u_ex       : callable (*coords_row) -> scalar (dim=1) or dim-tuple
    error_dict : ordered dict of {label: value} for error summary lines

    Column count adjusts with dim. For 1D the displacement column header is
    `ux_h` and the error header is `err_x` (uniform with 2D/3D).
    """
    coords = np.asarray(coords)
    if coords.ndim == 1:
        coords = coords.reshape(-1, 1)
    u_h = np.asarray(u_h)
    if u_h.ndim == 1:
        u_h = u_h.reshape(-1, 1)
    dim = coords.shape[1]
    assert u_h.shape[1] == dim, "u_h column count must match coord dim"
    axes = _AXES[:dim]

    coord_hdr = " | ".join(f"{ax:>10}"          for ax in axes)
    disp_hdr  = " | ".join(f"{'u'+ax+'_h':>15}" for ax in axes)
    exact_hdr = " | ".join(f"{'u'+ax+'_ex':>15}" for ax in axes)
    err_hdr   = " | ".join(f"{'err_'+ax:>15}"   for ax in axes)
    header    = f"{coord_hdr} | {disp_hdr} | {exact_hdr} | {err_hdr}"

    os.makedirs(results_dir, exist_ok=True)
    path = os.path.join(results_dir, f"{stem}.txt")
    with open(path, "w") as f:
        f.write(header + "\n")
        f.write("-" * len(header) + "\n")
        for c_row, uh_row in zip(coords, u_h):
            ue = u_ex(*c_row)
            ue = (ue,) if np.isscalar(ue) else tuple(ue)
            coord_str = " | ".join(f"{x:10.4f}"  for x in c_row)
            disp_str  = " | ".join(f"{x:15.6e}"  for x in uh_row)
            exact_str = " | ".join(f"{x:15.6e}"  for x in ue)
            err_str   = " | ".join(f"{abs(uh_row[k] - ue[k]):15.6e}"
                                   for k in range(dim))
            f.write(f"{coord_str} | {disp_str} | {exact_str} | {err_str}\n")
        f.write("\n")
        for label, val in error_dict.items():
            f.write(f"{label:12s} = {val:.6e}\n")


def write_convergence_table(stem, rows, results_dir):
    """
    Write convergence table to <results_dir>/<stem>.txt.

    rows : list of dicts with keys 'nx', 'h', and one key per error column.
           Rate columns are strings (empty for the first row).
    """
    os.makedirs(results_dir, exist_ok=True)
    path = os.path.join(results_dir, f"{stem}.txt")
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


def plot_convergence(stem, hs, series, title, results_dir, ylabel="Error"):
    """
    Save log-log convergence plot to <results_dir>/<stem>.png.

    series : list of {"label", "errors", "style"?} dicts
    Per-segment convergence rates are annotated above each line segment.
    """
    os.makedirs(results_dir, exist_ok=True)
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
    fig.savefig(os.path.join(results_dir, f"{stem}.png"), dpi=150)
    plt.close(fig)

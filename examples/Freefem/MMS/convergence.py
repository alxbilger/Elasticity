"""Shared convergence-study helpers for MMS 1D/2D/3D drivers."""

import numpy as np

from output import write_convergence_table


def run_convergence_series(nx_values, run_fn, h_fn, error_fns, banner,
                           results_dir, table_stem):
    """
    Run a mesh-refinement series and write a table.

    For each `nx` in `nx_values`, runs `run_fn(nx)` to get a solution,
    evaluates every error in `error_fns`, prints a live row, and writes
    the assembled table to ``<results_dir>/<table_stem>.txt``.

    Parameters
    ----------
    nx_values  : iterable of mesh-refinement parameters
    run_fn     : callable(nx) -> Solution (anything `error_fns` understands)
    h_fn       : callable(nx) -> mesh size h (drives the rate denominator
                 and the printed h column)
    error_fns  : ordered dict { label: callable(sol) -> float }
    banner     : string printed before the table header
    results_dir, table_stem : the file ``<results_dir>/<table_stem>.txt``
                 receives the assembled table

    Returns
    -------
    (hs, errors) : list of mesh sizes; dict {label: list of errors}.
    """
    hs     = []
    errors = {label: [] for label in error_fns}
    rows   = []
    err_labels = list(error_fns.keys())

    hdr = f"{'nx':>5} | {'h':>10}" + "".join(
        f" | {lbl:>14} | {('rate_' + lbl):>7}" for lbl in err_labels)
    print(f"\n{banner}\n{hdr}", flush=True)

    for k, nx in enumerate(nx_values):
        h   = h_fn(nx)
        sol = run_fn(nx)
        row = {"nx": nx, "h": h}
        for label, err_fn in error_fns.items():
            e_k  = err_fn(sol)
            rate = (f"{np.log(e_k / errors[label][-1]) / np.log(h / hs[-1]):.2f}"
                    if k > 0 else "")
            errors[label].append(e_k)
            row[label]           = e_k
            row[f"rate_{label}"] = rate

        line = f"{nx:5d} | {h:10.4f}"
        for lbl in err_labels:
            line += f" | {row[lbl]:14.6e} | {row[f'rate_{lbl}']:>7}"
        print(line, flush=True)

        hs.append(h)
        rows.append(row)

    write_convergence_table(table_stem, rows, results_dir)
    return hs, errors

"""Plot helpers shared across MMS 1D/2D drivers."""

import numpy as np


def annotate_convergence_rates(ax, hs, errors):
    """Annotate the per-segment log-log convergence rate above each segment."""
    hs     = np.asarray(hs)
    errors = np.asarray(errors)
    for k in range(1, len(hs)):
        rate  = np.log(errors[k] / errors[k-1]) / np.log(hs[k] / hs[k-1])
        x_mid = np.sqrt(hs[k] * hs[k-1])
        y_mid = np.sqrt(errors[k] * errors[k-1])
        ax.annotate(f"{rate:.2f}", xy=(x_mid, y_mid),
                    xytext=(0, 8), textcoords="offset points",
                    ha="center", fontsize=9)

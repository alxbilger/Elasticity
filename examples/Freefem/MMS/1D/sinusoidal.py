"""
Sinusoidal MMS (non-dimensional, x = x_dim/L ∈ [0,1], E = E_dim/L):

    u_ex(x) = sin(2πx)
    f(x)    = 4π²E·sin(2πx)

BC:
    u(0)   = 0                       (Dirichlet)
    u'(1)  = 2π·cos(2π) = 2π         (Neumann)  =>  F_N = E·u'(1) = 2πE
"""

import numpy as np

from manufactured_solution import MMSCase1D
from bar import line_quadrature, case_scene, run_reference_scene


class Sinusoidal(MMSCase1D):
    name              = "sinusoidal"
    plot_label        = r"$\sin(2\pi x)$"
    source_quadrature = staticmethod(line_quadrature(1))

    def u_ex(self, xi):
        return np.sin(2.0 * np.pi * xi)

    def du_ex(self, xi):
        return 2.0 * np.pi * np.cos(2.0 * np.pi * xi)

    def source(self, xi, E):
        return 4.0 * np.pi**2 * E * np.sin(2.0 * np.pi * xi)


mms = Sinusoidal()
createScene = case_scene(mms)


if __name__ == "__main__":
    run_reference_scene(mms)

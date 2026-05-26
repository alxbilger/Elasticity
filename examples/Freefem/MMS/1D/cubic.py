"""
Cubic MMS (non-dimensional, x = x_dim/L ∈ [0,1], E = E_dim/L):

    u_ex(x) = x²(1-x)
    f(x)    = E(6x - 2)

BC:
    u(0)   = 0     (Dirichlet)
    u'(1)  = -1    (Neumann)  =>  F_N = E·u'(1) = -E
"""

from manufactured_solution import MMSCase1D
from bar import line_quadrature, case_scene, run_reference_scene


class Cubic(MMSCase1D):
    name              = "cubic"
    plot_label        = r"$x^2(1-x)$"
    source_quadrature = staticmethod(line_quadrature(1))

    def u_ex(self, xi):
        return xi**2 * (1.0 - xi)

    def du_ex(self, xi):
        return 2.0 * xi - 3.0 * xi**2

    def source(self, xi, E):
        return E * (6.0 * xi - 2.0)


mms = Cubic()
createScene = case_scene(mms)


if __name__ == "__main__":
    run_reference_scene(mms)

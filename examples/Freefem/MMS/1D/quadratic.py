"""
Quadratic MMS (non-dimensional, x = x_dim/L ∈ [0,1], E = E_dim/L):

    u_ex(x) = x(1-x)
    f(x)    = 2E

BC:
    u(0)   = 0     (Dirichlet)
    u'(1)  = -1    (Neumann)  =>  F_N = E·u'(1) = -E
"""

from manufactured_solution import MMSCase1D
from bar import line_quadrature, case_scene, run_reference_scene


class Quadratic(MMSCase1D):
    name              = "quadratic"
    plot_label        = r"$x(1-x)$"
    source_quadrature = staticmethod(line_quadrature(1))

    def u_ex(self, xi):
        return xi * (1.0 - xi)

    def du_ex(self, xi):
        return 1.0 - 2.0 * xi

    def source(self, xi, E):
        return 2.0 * E


mms = Quadratic()
createScene = case_scene(mms)


if __name__ == "__main__":
    run_reference_scene(mms)

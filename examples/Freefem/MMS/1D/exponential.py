"""
Exponential MMS (non-dimensional, x = x_dim/L ∈ [0,1], E = E_dim/L):
    
    u_ex(x) = exp(x) - 1              du/dx   = exp(x)

    f(x)    = -E * exp(x)     
BC:
    u(0)   = 0     (Dirichlet)
    u'(1)  = e     (Neumann)  =>  F_N = E·u'(1) = E·e
"""
from manufactured_solution import MMSCase1D
from bar import line_quadrature, case_scene, run_reference_scene
import numpy as np

class Exponential(MMSCase1D):
    name              = "exponential"
    plot_label        = r"$e^x - 1$"
    source_quadrature = staticmethod(line_quadrature(1))

    def u_ex(self, xi):
        return np.exp(xi) - 1.0

    def du_ex(self, xi):
        return np.exp(xi)

    def source(self, xi, E):
        return -E * np.exp(xi)   # f = -E·u'' = -E·exp(x)

mms = Exponential()
createScene = case_scene(mms)

if __name__ == "__main__":
    run_reference_scene(mms)

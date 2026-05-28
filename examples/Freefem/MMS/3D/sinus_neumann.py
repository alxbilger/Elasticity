"""
Sinusoidal 3D MMS on [0,L]^3 with linear-elasticity constitutive law:

    u_ex(x, y, z) = A * ( sin(pi x / L) sin(pi y / L),
                          sin(pi y / L) sin(pi z / L),
                          sin(pi z / L) sin(pi x / L) )

    sigma = lambda tr(eps) I + 2 mu eps (3D / full Hooke).

BCs: minimal Dirichlet that removes the six rigid-body modes —
    - corner (0,0,0): pin all three components.
    - corner (L,0,0): pin u_y and u_z (kill rotation about z and y).
    - corner (0,L,0): pin u_z         (kill rotation about x).
Everything else is Neumann on the six faces (assembled by the framework
from `self.traction`).
"""

import numpy as np

from manufactured_solution import MMSCase3D, lame
from solid import (case_scene, run_reference_scene,
                   element_hex, hex_q1_rule)


SINUS_AMPLITUDE = 1e-1


class SinusNeumann(MMSCase3D):
    name       = "sinus_neumann"
    plot_label = (r"$u = A(\sin(\pi x/L)\sin(\pi y/L),\ "
                  r"\sin(\pi y/L)\sin(\pi z/L),\ "
                  r"\sin(\pi z/L)\sin(\pi x/L))$")

    source_quadrature_hex = staticmethod(hex_q1_rule(2))

    def u_ex(self, x, y, z, L):
        k = np.pi / L
        A = SINUS_AMPLITUDE
        return (A * np.sin(k*x) * np.sin(k*y),
                A * np.sin(k*y) * np.sin(k*z),
                A * np.sin(k*z) * np.sin(k*x))

    def grad_u_ex(self, x, y, z, L):
        k = np.pi / L
        A = SINUS_AMPLITUDE
        zero = 0.0 if np.isscalar(x) else np.zeros_like(np.asarray(x, float))
        dux_dx =  A * k * np.cos(k*x) * np.sin(k*y)
        dux_dy =  A * k * np.sin(k*x) * np.cos(k*y)
        dux_dz =  zero
        duy_dx =  zero
        duy_dy =  A * k * np.cos(k*y) * np.sin(k*z)
        duy_dz =  A * k * np.sin(k*y) * np.cos(k*z)
        duz_dx =  A * k * np.cos(k*x) * np.sin(k*z)
        duz_dy =  zero
        duz_dz =  A * k * np.cos(k*z) * np.sin(k*x)
        return np.array([[dux_dx, dux_dy, dux_dz],
                         [duy_dx, duy_dy, duy_dz],
                         [duz_dx, duz_dy, duz_dz]])

    def source(self, x, y, z, E, nu, L):
        lam, mu = lame(E, nu)
        A = SINUS_AMPLITUDE
        p = np.pi / L
        sx, sy, sz = np.sin(p*x), np.sin(p*y), np.sin(p*z)
        cx, cy, cz = np.cos(p*x), np.cos(p*y), np.cos(p*z)

        # Component-wise Laplacian
        lap_ux = -2 * p**2 * sx * sy
        lap_uy = -2 * p**2 * sy * sz
        lap_uz = -2 * p**2 * sz * sx

        # Gradient of the divergence
        d_divu_dx = p**2 * (-sx*sy + cz*cx)
        d_divu_dy = p**2 * ( cx*cy - sy*sz)
        d_divu_dz = p**2 * ( cy*cz - sz*sx)

        # Cauchy momentum: -div sigma = f, with sigma = lam tr(eps) I + 2 mu eps
        # for incompressible-style decomposition -((lam+mu) grad(div u) + mu lap u)
        fx = A * (-(lam + mu) * d_divu_dx - mu * lap_ux)
        fy = A * (-(lam + mu) * d_divu_dy - mu * lap_uy)
        fz = A * (-(lam + mu) * d_divu_dz - mu * lap_uz)
        return (fx, fy, fz)

    def apply_bcs(self, Solid, nodes_3d, L):
        eps = 1e-10
        xyz = nodes_3d[:, :3]

        def find_corner(pred):
            for k, (xk, yk, zk) in enumerate(xyz):
                if pred(xk, yk, zk):
                    return k
            raise RuntimeError("sinus_neumann: BC corner not found")

        i_origin = find_corner(lambda x, y, z: x < eps     and y < eps     and z < eps)
        i_xL     = find_corner(lambda x, y, z: x > L - eps and y < eps     and z < eps)
        i_yL     = find_corner(lambda x, y, z: x < eps     and y > L - eps and z < eps)

        Solid.addObject("FixedProjectiveConstraint",
                        name="fix_origin", indices=i_origin)
        Solid.addObject("PartialFixedProjectiveConstraint",
                        name="fix_x_axis", template="Vec3d",
                        indices=i_xL, fixedDirections="0 1 1")
        Solid.addObject("PartialFixedProjectiveConstraint",
                        name="fix_y_axis", template="Vec3d",
                        indices=i_yL, fixedDirections="0 0 1")


mms         = SinusNeumann()
createScene = case_scene(mms, element_hex)


if __name__ == "__main__":
    run_reference_scene(element_hex, mms)

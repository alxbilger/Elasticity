"""
Incompressible 2D MMS on [0,L]^2 with linear-elasticity constitutive law:

    u_ex(x, y) = ( sin(pi x / L) cos(pi y / L),
                  -cos(pi x / L) sin(pi y / L) )

    div u = 0 everywhere — the lambda * tr(eps) term in sigma vanishes
    pointwise, so the discretization is exercised in the lambda → infinity
    (nu → 0.5) limit. Useful for probing volumetric locking.

    BCs: minimal Dirichlet (corner anchor + rotation-killer) so the system is
    constrained only by the Neumann tractions on the four edges — that's the
    point of the test.
"""

import numpy as np

from manufactured_solution import MMSCase2D, lame
from beam import (case_scene, run_reference_scene,
                  element_quad, element_tri,
                  quad_q1_rule, tri_p1_rule)


class Incompressible(MMSCase2D):
    name       = "incompressible"
    plot_label = (r"$u_x = \sin(\pi x/L)\cos(\pi y/L),\ "
                  r"u_y = -\cos(\pi x/L)\sin(\pi y/L)$")

    source_quadrature_quad = staticmethod(quad_q1_rule(2))
    source_quadrature_tri  = staticmethod(tri_p1_rule(3))

    def u_ex(self, x, y, L):
        k = np.pi / L
        return ( np.sin(k * x) * np.cos(k * y),
                -np.cos(k * x) * np.sin(k * y))

    def grad_u_ex(self, x, y, L):
        k = np.pi / L
        dux_dx =  k * np.cos(k * x) * np.cos(k * y)
        dux_dy = -k * np.sin(k * x) * np.sin(k * y)
        duy_dx =  k * np.sin(k * x) * np.sin(k * y)
        duy_dy = -k * np.cos(k * x) * np.cos(k * y)
        return np.array([[dux_dx, dux_dy],
                         [duy_dx, duy_dy]])

    def source(self, x, y, E, nu, L, dim):
        lam, mu = lame(E, nu, dim)
        k  = np.pi / L
        ux =  np.sin(k * x) * np.cos(k * y)
        uy = -np.cos(k * x) * np.sin(k * y)
        d2ux_dxx = -k**2 * ux
        d2ux_dyy = -k**2 * ux
        d2ux_dxy = -k**2 * np.cos(k * x) * np.sin(k * y)
        d2uy_dxx = -k**2 * uy
        d2uy_dyy = -k**2 * uy
        d2uy_dxy =  k**2 * np.sin(k * x) * np.cos(k * y)
        fx = -((lam + 2*mu) * d2ux_dxx + lam * d2uy_dxy
             + mu * (d2ux_dyy + d2uy_dxy))
        fy = -(mu * (d2ux_dxy + d2uy_dxx) + lam * d2ux_dxy
             + (lam + 2*mu) * d2uy_dyy)
        return (fx, fy)

    def apply_bcs(self, Beam, nodes_2d, L, dim):
        # Minimal Dirichlet: pin (0,0) in both directions to kill translation,
        # pin uy at (L,0) to kill rotation. Everything else is Neumann.
        tmpl = "Vec3d" if dim == "3d" else "Vec2d"
        xy   = nodes_2d[:, :2]
        eps  = 1e-10
        idx_corner = [k for k, (xk, yk) in enumerate(xy)
                      if xk < eps and yk < eps]
        idx_rot    = [k for k, (xk, yk) in enumerate(xy)
                      if xk > L - eps and yk < eps]

        if dim == "2d":
            Beam.addObject("PartialFixedProjectiveConstraint",
                           name="fix_corner", template=tmpl,
                           indices=" ".join(map(str, idx_corner)),
                           fixedDirections="1 1")
            Beam.addObject("PartialFixedProjectiveConstraint",
                           name="fix_rot", template=tmpl,
                           indices=" ".join(map(str, idx_rot)),
                           fixedDirections="0 1")
        else:
            all_idx = " ".join(map(str, range(len(nodes_2d))))
            Beam.addObject("PartialFixedProjectiveConstraint",
                           name="fix_corner", template=tmpl,
                           indices=" ".join(map(str, idx_corner)),
                           fixedDirections="1 1 0")
            Beam.addObject("PartialFixedProjectiveConstraint",
                           name="fix_rot", template=tmpl,
                           indices=" ".join(map(str, idx_rot)),
                           fixedDirections="0 1 0")
            Beam.addObject("PartialFixedProjectiveConstraint",
                           name="fix_z", template=tmpl,
                           indices=all_idx,
                           fixedDirections="0 0 1")


mms         = Incompressible()
createScene = case_scene(mms, element_quad)


if __name__ == "__main__":
    run_reference_scene(element_quad, mms)

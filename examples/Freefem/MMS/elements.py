"""Element strategies for MMS scenes.

Each strategy declares:

    LABEL              : human-readable identifier (used in plot legends)
    ELEMENT_RULE       : an element rule from `fem.py` driving the L²/H¹ norms
    _source_rule(mms)  : the body-force quadrature rule (read off the mms;
                         raises if the case did not set it)
    add_topology(parent_node)   : add the SOFA topology container and return it
    read_connectivity(topology) : extract the connectivity array post-init
    compute_nodal_forces(...)   : assemble body + traction nodal forces
    compute_l2(sol, mms, L)     : L² error norm on a solution
    compute_h1(sol, mms, L)     : H¹ semi-norm error on a solution

Two `_ElementBase` classes live here — one for 2D (4 boundary edges,
dim-dependent plane stress / plane strain) and one for 3D (6 boundary
quads, single constitutive branch). Concrete elements (`_QuadElement`,
`_TriElement`, `_HexElement`) subclass the matching base and pin the
topology container choice and the rules.
"""

import numpy as np

from fem import (
    assemble_nodal_forces,
    assemble_traction,
    l2_error,
    h1_semi_error,
    quad_q1_rule,
    tri_p1_rule,
    hex_q1_rule,
    edge_line_rule,
    quad_face_rule,
)


# ---------------------------------------------------------------------------
# Boundary-facet helpers
# ---------------------------------------------------------------------------

def _boundary_edges(nx, ny):
    """Return (bottom, top, left, right) edge lists for a structured nx×ny grid."""
    bottom = [(i, i + 1)                                    for i in range(nx - 1)]
    top    = [((ny - 1) * nx + i, (ny - 1) * nx + i + 1)    for i in range(nx - 1)]
    left   = [(j * nx, (j + 1) * nx)                        for j in range(ny - 1)]
    right  = [(j * nx + (nx - 1), (j + 1) * nx + (nx - 1))  for j in range(ny - 1)]
    return bottom, top, left, right


def _boundary_quads(nx, ny, nz):
    """Return (xm, xp, ym, yp, zm, zp) quad lists for a structured nx×ny×nz grid.

    Each quad is a 4-tuple of node indices using SOFA's regular-grid index
    convention `idx(i, j, k) = i + j*nx + k*nx*ny`. Orientation is consistent
    per face but its sign does not affect `assemble_traction`, which uses
    |t_xi × t_eta| for the surface area and takes the outward normal from
    the caller.
    """
    def idx(i, j, k):
        return i + j * nx + k * nx * ny

    xm = [(idx(0, j, k),    idx(0, j+1, k),    idx(0, j+1, k+1),    idx(0, j, k+1))
          for k in range(nz - 1) for j in range(ny - 1)]
    xp = [(idx(nx-1, j, k), idx(nx-1, j+1, k), idx(nx-1, j+1, k+1), idx(nx-1, j, k+1))
          for k in range(nz - 1) for j in range(ny - 1)]
    ym = [(idx(i, 0, k),    idx(i+1, 0, k),    idx(i+1, 0, k+1),    idx(i, 0, k+1))
          for k in range(nz - 1) for i in range(nx - 1)]
    yp = [(idx(i, ny-1, k), idx(i+1, ny-1, k), idx(i+1, ny-1, k+1), idx(i, ny-1, k+1))
          for k in range(nz - 1) for i in range(nx - 1)]
    zm = [(idx(i, j, 0),    idx(i+1, j, 0),    idx(i+1, j+1, 0),    idx(i, j+1, 0))
          for j in range(ny - 1) for i in range(nx - 1)]
    zp = [(idx(i, j, nz-1), idx(i+1, j, nz-1), idx(i+1, j+1, nz-1), idx(i, j+1, nz-1))
          for j in range(ny - 1) for i in range(nx - 1)]
    return xm, xp, ym, yp, zm, zp


# ---------------------------------------------------------------------------
# 2D elements (Vec2d / plane stress  OR  Vec3d / plane strain via dim arg)
# ---------------------------------------------------------------------------

class _ElementBase2D:
    @classmethod
    def compute_nodal_forces(cls, nodes_2d, conn, mms, L, E, nu, nx, ny, dim):
        xy = nodes_2d[:, :2]

        F = assemble_nodal_forces(
            lambda x, y: mms.source(x, y, E, nu, L, dim),
            xy, conn, cls._source_rule(mms))

        bottom, top, left, right = _boundary_edges(nx, ny)
        sides = [(bottom, 0.0, -1.0),
                 (top,    0.0, +1.0),
                 (left,  -1.0,  0.0),
                 (right, +1.0,  0.0)]
        edge_rule = edge_line_rule(2)
        for edges, nrm_x, nrm_y in sides:
            F += assemble_traction(
                lambda x, y, nx=nrm_x, ny=nrm_y:
                    mms.traction(x, y, nx, ny, E, nu, L, dim),
                xy, edges, edge_rule)
        return F

    @classmethod
    def compute_l2(cls, sol, mms, L):
        xy = sol.nodes[:, :2]
        return l2_error(
            xy, sol.conn, np.column_stack([sol.ux, sol.uy]),
            lambda x, y: mms.u_ex(x, y, L),
            cls.ELEMENT_RULE)

    @classmethod
    def compute_h1(cls, sol, mms, L):
        xy = sol.nodes[:, :2]
        return h1_semi_error(
            xy, sol.conn, np.column_stack([sol.ux, sol.uy]),
            lambda x, y: mms.grad_u_ex(x, y, L),
            cls.ELEMENT_RULE)


class _QuadElement(_ElementBase2D):
    LABEL        = "Q1 quad"
    ELEMENT_RULE = staticmethod(quad_q1_rule(2))   # used for L²/H¹ error norms

    @staticmethod
    def _source_rule(mms):
        rule = mms.source_quadrature_quad
        if rule is None:
            raise ValueError(
                f"{type(mms).__name__}.source_quadrature_quad must be set")
        return rule

    @staticmethod
    def add_topology(Beam):
        topology = Beam.addObject("QuadSetTopologyContainer", name="topology",
                                  quads="@../Grid/grid.quads",
                                  position="@../Grid/grid.position")
        Beam.addObject("QuadSetTopologyModifier")
        return topology

    @staticmethod
    def read_connectivity(topology):
        return topology.quads.array().copy()

    @staticmethod
    def to_triangles(conn):
        tris = []
        for q in conn:
            tris.append([q[0], q[1], q[2]])
            tris.append([q[0], q[2], q[3]])
        return np.array(tris)


class _TriElement(_ElementBase2D):
    LABEL        = "P1 tri"
    ELEMENT_RULE = staticmethod(tri_p1_rule(3))    # used for L²/H¹ error norms

    @staticmethod
    def _source_rule(mms):
        rule = mms.source_quadrature_tri
        if rule is None:
            raise ValueError(
                f"{type(mms).__name__}.source_quadrature_tri must be set")
        return rule

    @staticmethod
    def add_topology(Beam):
        topology = Beam.addObject("TriangleSetTopologyContainer", name="topology")
        Beam.addObject("Quad2TriangleTopologicalMapping",
                       input="@../Grid/grid", output="@topology")
        Beam.addObject("TriangleSetTopologyModifier")
        return topology

    @staticmethod
    def read_connectivity(topology):
        return topology.triangles.array().copy()

    @staticmethod
    def to_triangles(conn):
        return conn


# ---------------------------------------------------------------------------
# 3D elements (Vec3d / full Hooke)
# ---------------------------------------------------------------------------

class _ElementBase3D:
    @classmethod
    def compute_nodal_forces(cls, nodes_3d, conn, mms, L, E, nu, nx, ny, nz):
        xyz = nodes_3d[:, :3]

        F = assemble_nodal_forces(
            lambda x, y, z: mms.source(x, y, z, E, nu, L),
            xyz, conn, cls._source_rule(mms))

        xm, xp, ym, yp, zm, zp = _boundary_quads(nx, ny, nz)
        sides = [(xm, -1.0, 0.0, 0.0),
                 (xp, +1.0, 0.0, 0.0),
                 (ym,  0.0, -1.0, 0.0),
                 (yp,  0.0, +1.0, 0.0),
                 (zm,  0.0, 0.0, -1.0),
                 (zp,  0.0, 0.0, +1.0)]
        face_rule = quad_face_rule(2)
        for quads, nrm_x, nrm_y, nrm_z in sides:
            F += assemble_traction(
                lambda x, y, z, nx=nrm_x, ny=nrm_y, nz=nrm_z:
                    mms.traction(x, y, z, nx, ny, nz, E, nu, L),
                xyz, quads, face_rule)
        return F

    @classmethod
    def compute_l2(cls, sol, mms, L):
        return l2_error(
            sol.nodes, sol.conn,
            np.column_stack([sol.ux, sol.uy, sol.uz]),
            lambda x, y, z: mms.u_ex(x, y, z, L),
            cls.ELEMENT_RULE)

    @classmethod
    def compute_h1(cls, sol, mms, L):
        return h1_semi_error(
            sol.nodes, sol.conn,
            np.column_stack([sol.ux, sol.uy, sol.uz]),
            lambda x, y, z: mms.grad_u_ex(x, y, z, L),
            cls.ELEMENT_RULE)


class _HexElement(_ElementBase3D):
    LABEL        = "Q1 hex"
    ELEMENT_RULE = staticmethod(hex_q1_rule(2))   # used for L²/H¹ error norms

    @staticmethod
    def _source_rule(mms):
        rule = mms.source_quadrature_hex
        if rule is None:
            raise ValueError(
                f"{type(mms).__name__}.source_quadrature_hex must be set")
        return rule

    @staticmethod
    def add_topology(Solid):
        topology = Solid.addObject("HexahedronSetTopologyContainer",
                                   name="topology",
                                   hexahedra="@../Grid/grid.hexahedra",
                                   position="@../Grid/grid.position")
        Solid.addObject("HexahedronSetTopologyModifier")
        return topology

    @staticmethod
    def read_connectivity(topology):
        return topology.hexahedra.array().copy()


# ---------------------------------------------------------------------------
# Instances (one per element type)
# ---------------------------------------------------------------------------

element_quad = _QuadElement()
element_tri  = _TriElement()
element_hex  = _HexElement()

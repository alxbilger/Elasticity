import json
import os
import numpy as np
import matplotlib.pyplot as plt
import Sofa
import Sofa.Core
import Sofa.Simulation

RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results_3d")
SINUS_AMPLITUDE = 1e-1
 

def load_params(path=None):
    if path is None:
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "params.json")
    with open(path) as f:
        return json.load(f)


def lame(E, nu):
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu  = E / (2 * (1 + nu))
    return lam, mu


def node_idx(i, j, k, nx, ny):
    return i + j * nx + k * nx * ny


def _build_element_connectivity(nx, ny, nz):
    i = np.arange(nx-1); j = np.arange(ny-1); k = np.arange(nz-1)
    ii, jj, kk = np.meshgrid(i, j, k, indexing='ij')
    ii, jj, kk = ii.ravel(), jj.ravel(), kk.ravel()
    return np.column_stack([
        ii   + jj*nx     + kk*nx*ny,
        ii+1 + jj*nx     + kk*nx*ny,
        ii+1 + (jj+1)*nx + kk*nx*ny,
        ii   + (jj+1)*nx + kk*nx*ny,
        ii   + jj*nx     + (kk+1)*nx*ny,
        ii+1 + jj*nx     + (kk+1)*nx*ny,
        ii+1 + (jj+1)*nx + (kk+1)*nx*ny,
        ii   + (jj+1)*nx + (kk+1)*nx*ny,
    ])  # (n_e, 8)

 
def u_ex_linear(x, y, z, L):
    return x / L, y / L, z / L


def stress_linear(E, nu, L):
    lam, mu = lame(E, nu)
    eps = np.diag([1/L, 1/L, 1/L])
    return lam * np.trace(eps) * np.eye(3) + 2 * mu * eps


def traction_linear(x, y, z, nx_c, ny_c, nz_c, E, nu, L):
    t = stress_linear(E, nu, L) @ np.array([nx_c, ny_c, nz_c])
    return t[0], t[1], t[2]

 

def u_ex_sinus_neumann(x, y, z, L):
    return (SINUS_AMPLITUDE * np.sin(np.pi*x/L)*np.sin(np.pi*y/L),
            SINUS_AMPLITUDE * np.sin(np.pi*y/L)*np.sin(np.pi*z/L),
            SINUS_AMPLITUDE * np.sin(np.pi*z/L)*np.sin(np.pi*x/L))

def _sn_dux_dx(x,y,z,L): return (np.pi/L)*np.cos(np.pi*x/L)*np.sin(np.pi*y/L)
def _sn_dux_dy(x,y,z,L): return (np.pi/L)*np.sin(np.pi*x/L)*np.cos(np.pi*y/L)
def _sn_dux_dz(x,y,z,L): return np.zeros_like(x) if hasattr(x,'__len__') else 0.0
def _sn_duy_dx(x,y,z,L): return np.zeros_like(x) if hasattr(x,'__len__') else 0.0
def _sn_duy_dy(x,y,z,L): return (np.pi/L)*np.cos(np.pi*y/L)*np.sin(np.pi*z/L)
def _sn_duy_dz(x,y,z,L): return (np.pi/L)*np.sin(np.pi*y/L)*np.cos(np.pi*z/L)
def _sn_duz_dx(x,y,z,L): return (np.pi/L)*np.sin(np.pi*z/L)*np.cos(np.pi*x/L)
def _sn_duz_dy(x,y,z,L): return np.zeros_like(x) if hasattr(x,'__len__') else 0.0
def _sn_duz_dz(x,y,z,L): return (np.pi/L)*np.cos(np.pi*z/L)*np.sin(np.pi*x/L)

def stress_sinus(x, y, z, E, nu, L):
    lam, mu = lame(E, nu)
    exx = _sn_dux_dx(x,y,z,L); eyy = _sn_duy_dy(x,y,z,L); ezz = _sn_duz_dz(x,y,z,L)
    exy = 0.5*(_sn_dux_dy(x,y,z,L) + _sn_duy_dx(x,y,z,L))
    exz = 0.5*(_sn_dux_dz(x,y,z,L) + _sn_duz_dx(x,y,z,L))
    eyz = 0.5*(_sn_duy_dz(x,y,z,L) + _sn_duz_dy(x,y,z,L))
    tr  = exx + eyy + ezz
    return SINUS_AMPLITUDE * np.array([
        [lam*tr+2*mu*exx, 2*mu*exy,        2*mu*exz       ],
        [2*mu*exy,        lam*tr+2*mu*eyy, 2*mu*eyz       ],
        [2*mu*exz,        2*mu*eyz,        lam*tr+2*mu*ezz],
    ])

def traction_sinus(x, y, z, nx_c, ny_c, nz_c, E, nu, L):
    t = stress_sinus(x, y, z, E, nu, L) @ np.array([nx_c, ny_c, nz_c])
    return t[0], t[1], t[2]

def body_force_sinus(x, y, z, E, nu, L):
    lam, mu = lame(E, nu)
    p = np.pi / L
    sx,sy,sz = np.sin(p*x), np.sin(p*y), np.sin(p*z)
    cx,cy,cz = np.cos(p*x), np.cos(p*y), np.cos(p*z)
    lap_ux = -2*p**2 * sx*sy
    lap_uy = -2*p**2 * sy*sz
    lap_uz = -2*p**2 * sz*sx
    d_divu_dx = p**2 * (-sx*sy + cz*cx)
    d_divu_dy = p**2 * ( cx*cy - sy*sz)
    d_divu_dz = p**2 * ( cy*cz - sz*sx)
    fx = SINUS_AMPLITUDE * (-(lam+mu)*d_divu_dx - mu*lap_ux)
    fy = SINUS_AMPLITUDE * (-(lam+mu)*d_divu_dy - mu*lap_uy)
    fz = SINUS_AMPLITUDE * (-(lam+mu)*d_divu_dz - mu*lap_uz)
    return fx, fy, fz

 

def u_ex_sinus_dirichlet(x, y, z, L):
    val = SINUS_AMPLITUDE * np.sin(np.pi*x/L)*np.sin(np.pi*y/L)*np.sin(np.pi*z/L)
    return val, val, val

def _sd_dx(x,y,z,L): return (np.pi/L)*np.cos(np.pi*x/L)*np.sin(np.pi*y/L)*np.sin(np.pi*z/L)
def _sd_dy(x,y,z,L): return (np.pi/L)*np.sin(np.pi*x/L)*np.cos(np.pi*y/L)*np.sin(np.pi*z/L)
def _sd_dz(x,y,z,L): return (np.pi/L)*np.sin(np.pi*x/L)*np.sin(np.pi*y/L)*np.cos(np.pi*z/L)

def body_force_sinus_dirichlet(x, y, z, E, nu, L):
    lam, mu = lame(E, nu)
    p = np.pi / L
    sx,sy,sz = np.sin(p*x), np.sin(p*y), np.sin(p*z)
    cx,cy,cz = np.cos(p*x), np.cos(p*y), np.cos(p*z)

    lap_s3 = -3*p**2 * sx*sy*sz
    d_divu_dx = p**2*(-sx*sy*sz + cx*cy*sz + cx*sy*cz)
    d_divu_dy = p**2*( cx*cy*sz - sx*sy*sz + sx*cy*cz)
    d_divu_dz = p**2*( cx*sy*cz + sx*cy*cz - sx*sy*sz)
    fx = SINUS_AMPLITUDE * (-(lam+mu)*d_divu_dx - mu*lap_s3)
    fy = SINUS_AMPLITUDE * (-(lam+mu)*d_divu_dy - mu*lap_s3)
    fz = SINUS_AMPLITUDE * (-(lam+mu)*d_divu_dz - mu*lap_s3)
    return fx, fy, fz

def get_u_ex(mode):
    return {'linear_neumann':  u_ex_linear,
            'sinus_neumann':   u_ex_sinus_neumann,
            'sinus_dirichlet': u_ex_sinus_dirichlet}[mode]

def get_derivatives(mode):
    if mode == 'linear_neumann':
        one = lambda x,y,z,L: np.ones_like(x)/L if hasattr(x,'__len__') else 1/L
        zer = lambda x,y,z,L: np.zeros_like(x)  if hasattr(x,'__len__') else 0.0
        return {'dux_dx':one,'dux_dy':zer,'dux_dz':zer,
                'duy_dx':zer,'duy_dy':one,'duy_dz':zer,
                'duz_dx':zer,'duz_dy':zer,'duz_dz':one}
    elif mode == 'sinus_neumann':
        def _s(fn): return lambda x,y,z,L,_f=fn: _f(x,y,z,L)*SINUS_AMPLITUDE
        return {'dux_dx':_s(_sn_dux_dx),'dux_dy':_s(_sn_dux_dy),'dux_dz':_s(_sn_dux_dz),
                'duy_dx':_s(_sn_duy_dx),'duy_dy':_s(_sn_duy_dy),'duy_dz':_s(_sn_duy_dz),
                'duz_dx':_s(_sn_duz_dx),'duz_dy':_s(_sn_duz_dy),'duz_dz':_s(_sn_duz_dz)}
    else:  # sinus_dirichlet : ux=uy=uz=sin*sin*sin
        def _s(fn): return lambda x,y,z,L,_f=fn: _f(x,y,z,L)*SINUS_AMPLITUDE
        sd = _s(_sd_dx); se = _s(_sd_dy); sf = _s(_sd_dz)
        return {'dux_dx':sd,'dux_dy':se,'dux_dz':sf,
                'duy_dx':sd,'duy_dy':se,'duy_dz':sf,
                'duz_dx':sd,'duz_dy':se,'duz_dz':sf}

 

def _hexa_shape(xi, eta, zeta):
    N = 0.125 * np.array([
        (1-xi)*(1-eta)*(1-zeta), (1+xi)*(1-eta)*(1-zeta),
        (1+xi)*(1+eta)*(1-zeta), (1-xi)*(1+eta)*(1-zeta),
        (1-xi)*(1-eta)*(1+zeta), (1+xi)*(1-eta)*(1+zeta),
        (1+xi)*(1+eta)*(1+zeta), (1-xi)*(1+eta)*(1+zeta),
    ])
    dN_dxi = 0.125 * np.array([
        -(1-eta)*(1-zeta), (1-eta)*(1-zeta), (1+eta)*(1-zeta),-(1+eta)*(1-zeta),
        -(1-eta)*(1+zeta), (1-eta)*(1+zeta), (1+eta)*(1+zeta),-(1+eta)*(1+zeta),
    ])
    dN_deta = 0.125 * np.array([
        -(1-xi)*(1-zeta),-(1+xi)*(1-zeta), (1+xi)*(1-zeta), (1-xi)*(1-zeta),
        -(1-xi)*(1+zeta),-(1+xi)*(1+zeta), (1+xi)*(1+zeta), (1-xi)*(1+zeta),
    ])
    dN_dzeta = 0.125 * np.array([
        -(1-xi)*(1-eta),-(1+xi)*(1-eta),-(1+xi)*(1+eta),-(1-xi)*(1+eta),
         (1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta),
    ])
    return N, dN_dxi, dN_deta, dN_dzeta


# Precomputed shape functions and derivatives at all 8 Gauss points (2-point rule, gw=1)
_GP = np.array([-1/np.sqrt(3), 1/np.sqrt(3)])
_N_GP  = np.empty((8, 8))    # (n_gp, n_nodes)
_dN_GP = np.empty((8, 3, 8)) # (n_gp, 3, n_nodes)
for _g, (_xi, _eta, _zeta) in enumerate(
        [(_xi, _eta, _zeta) for _xi in _GP for _eta in _GP for _zeta in _GP]):
    _N, _dNxi, _dNeta, _dNzeta = _hexa_shape(_xi, _eta, _zeta)
    _N_GP[_g]  = _N
    _dN_GP[_g] = np.array([_dNxi, _dNeta, _dNzeta])


def _batch_det3(J):
    """Analytical determinant for a batch of (n, 3, 3) matrices."""
    return (J[:,0,0] * (J[:,1,1]*J[:,2,2] - J[:,1,2]*J[:,2,1])
          - J[:,0,1] * (J[:,1,0]*J[:,2,2] - J[:,1,2]*J[:,2,0])
          + J[:,0,2] * (J[:,1,0]*J[:,2,1] - J[:,1,1]*J[:,2,0]))


def _batch_detinv3(J):
    """Analytical det and inverse for (n, 3, 3) matrices via cofactors."""
    a=J[:,0,0]; b=J[:,0,1]; c=J[:,0,2]
    d=J[:,1,0]; e=J[:,1,1]; f=J[:,1,2]
    g=J[:,2,0]; h=J[:,2,1]; k=J[:,2,2]
    c00=e*k-f*h;  c10=f*g-d*k;  c20=d*h-e*g
    c01=c*h-b*k;  c11=a*k-c*g;  c21=b*g-a*h
    c02=b*f-c*e;  c12=c*d-a*f;  c22=a*e-b*d
    det = a*c00 + b*c10 + c*c20
    s = 1.0 / det
    Jinv = np.empty_like(J)
    Jinv[:,0,0]=c00*s; Jinv[:,0,1]=c01*s; Jinv[:,0,2]=c02*s
    Jinv[:,1,0]=c10*s; Jinv[:,1,1]=c11*s; Jinv[:,1,2]=c12*s
    Jinv[:,2,0]=c20*s; Jinv[:,2,1]=c21*s; Jinv[:,2,2]=c22*s
    return det, Jinv


def compute_body_force_nodal(nodes, nx, ny, nz, E, nu, L, body_force_fn):
    all_idx = _build_element_connectivity(nx, ny, nz)  # (n_e, 8)
    xe_all  = nodes[all_idx]                            # (n_e, 8, 3)
    F = np.zeros((len(nodes), 3))

    for g in range(8):
        J_all    = _dN_GP[g] @ xe_all                            # (n_e, 3, 3)
        detJ_all = _batch_det3(J_all)                            # (n_e,)
        xg_all   = xe_all.transpose(0, 2, 1) @ _N_GP[g]         # (n_e, 3)

        fx, fy, fz = body_force_fn(xg_all[:,0], xg_all[:,1], xg_all[:,2], E, nu, L)
        f_all  = np.stack([fx, fy, fz], axis=1)                  # (n_e, 3)
        w_all  = np.abs(detJ_all)                                 # gw=1 for all points

        # contrib[e, local_node, d] = N_g[local_node] * f_all[e,d] * w_all[e]
        contrib = (_N_GP[g][np.newaxis, :, np.newaxis]
                   * f_all[:, np.newaxis, :]
                   * w_all[:, np.newaxis, np.newaxis])            # (n_e, 8, 3)
        for d in range(3):
            F[:, d] += np.bincount(all_idx.ravel(),
                                   weights=contrib[:, :, d].ravel(),
                                   minlength=len(nodes))
    return F

 

def compute_surface_traction_nodal(nodes, quads, E, nu, L, traction_fn):
    
    gp = np.array([-1/np.sqrt(3), 1/np.sqrt(3)])
    gw = np.array([1.0, 1.0])
    F  = np.zeros((len(nodes), 3))
    centroid = nodes.mean(axis=0)

    for quad in quads:
        x = nodes[quad]  
        dN_dxi0  = np.array([-1.,  1.,  1., -1.]) / 4
        dN_deta0 = np.array([-1., -1.,  1.,  1.]) / 4
        J0   = np.cross(dN_dxi0 @ x, dN_deta0 @ x)
        sign = 1.0 if np.dot(J0, x.mean(axis=0) - centroid) > 0 else -1.0

        for gi, xi in enumerate(gp):
            for gj, eta in enumerate(gp):
                N       = np.array([(1-xi)*(1-eta),(1+xi)*(1-eta),
                                    (1+xi)*(1+eta),(1-xi)*(1+eta)]) / 4
                dN_dxi  = np.array([-(1-eta), (1-eta), (1+eta),-(1+eta)]) / 4
                dN_deta = np.array([-(1-xi), -(1+xi),  (1+xi), (1-xi) ]) / 4
                J_vec   = np.cross(dN_dxi @ x, dN_deta @ x)
                Jmag    = np.linalg.norm(J_vec)
                n_hat   = sign * J_vec / Jmag
                xg      = N @ x
                T       = np.array(traction_fn(*xg, *n_hat, E, nu, L))
                F[quad] += np.outer(N, T) * (gw[gi]*gw[gj]*Jmag)

    residual = np.abs(F.sum(axis=0))
    return F
 

class MMSForceController(Sofa.Core.Controller):
    def __init__(self, dofs, cff, quad_topo, nx, ny, nz,
                 E, nu, L, mode, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.dofs      = dofs
        self.cff       = cff
        self.quad_topo = quad_topo
        self.nx, self.ny, self.nz = nx, ny, nz
        self.E, self.nu, self.L   = E, nu, L
        self.mode = mode

    def onSimulationInitDoneEvent(self, event):
        nodes = self.dofs.position.array().copy()
        quads = np.array(self.quad_topo.quads.array())
        E, nu, L       = self.E, self.nu, self.L
        nx, ny, nz     = self.nx, self.ny, self.nz
        mode           = self.mode

        if mode == 'linear_neumann':
            F = compute_surface_traction_nodal(
                    nodes, quads, E, nu, L, traction_linear)
        elif mode == 'sinus_neumann':
            F_body = compute_body_force_nodal(nodes, nx, ny, nz, E, nu, L, body_force_sinus)
            F_surf = compute_surface_traction_nodal(nodes, quads, E, nu, L, traction_sinus)
            print(f"  [body] résidu global F_body : {np.abs(F_body.sum(axis=0))}")
            print(f"  [total] résidu F_body+F_surf : {np.abs((F_body + F_surf).sum(axis=0))}")
            F = F_body + F_surf
        elif mode == 'sinus_dirichlet':
            F = compute_body_force_nodal(
                    nodes, nx, ny, nz, E, nu, L, body_force_sinus_dirichlet)
        else:
            raise ValueError(f"Unknown: {mode}")

        self.cff.forces.value = F.tolist()


# =========================== SOFA scene  ================================================ 

def build_scene_3d(rootNode, L=1.0, E=1e6, nu=0.3, nx=6, ny=6, nz=6,
                   mode='linear_neumann', visual=False):
    plugins = [
        "Elasticity",
        "Sofa.Component.Constraint.Projective",
        "Sofa.Component.Engine.Select",
        "Sofa.Component.LinearSolver.Direct",
        "Sofa.Component.MechanicalLoad",
        "Sofa.Component.ODESolver.Backward",
        "Sofa.Component.StateContainer",
        "Sofa.Component.Topology.Container.Dynamic",
        "Sofa.Component.Topology.Container.Grid",
        "Sofa.Component.Topology.Mapping",
        "Sofa.Component.LinearSolver.Iterative",
        "Sofa.Component.LinearSolver.Preconditioner",
        "Sofa.Component.LinearSystem",
        "Sofa.Component.Mapping.Linear"
    ]
    if visual:
        plugins += ["Sofa.Component.Visual","Sofa.GL.Component.Rendering3D"]

    rootNode.addObject('RequiredPlugin', pluginName=plugins)
    rootNode.addObject('DefaultAnimationLoop')
    if visual:
        rootNode.addObject('DefaultVisualManagerLoop')
        rootNode.addObject('VisualStyle',
                           displayFlags='showVisualModels showWireframe')
    rootNode.gravity.value = [0,0,0]
    rootNode.dt.value = 1.0

    regGrid = rootNode.addChild('GridTopology')
    regGrid.addObject('RegularGridTopology', name="grid",
                      nx=nx, ny=ny, nz=nz,
                      min=[0.,0.,0.], max=[L,L,L])

    solid = rootNode.addChild('Solid3D')
    solid.addObject('NewtonRaphsonSolver', name="newtonSolver",
                    printLog=False, warnWhenLineSearchFails=True
                    , maxNbIterationsNewton=2
                    , relativeSuccessiveStoppingThreshold=0
                    , absoluteResidualStoppingThreshold=0
                    , absoluteEstimateDifferenceThreshold=0
                    , relativeInitialStoppingThreshold=0
                    , relativeEstimateDifferenceThreshold=0)
    solid.addObject("MatrixLinearSystem", name="ssorSystem", template="CompressedRowSparseMatrixd")
    solid.addObject("SSORPreconditioner", name="preconditioner", printLog=False)
    solid.addObject("PreconditionedMatrixFreeSystem", name="solverSystem", template="GraphScattered", preconditionerSystem="@ssorSystem")
    solid.addObject("PCGLinearSolver", name="linearSolver", template="GraphScattered"
        , iterations=3000, tolerance=1e-15, preconditioner="@preconditioner", printLog=True)
    #solid.addObject('SparseLDLSolver', name="linearSolver", template="CompressedRowSparseMatrixd")
    solid.addObject('StaticSolver', name="staticSolver",
                    newtonSolver="@newtonSolver", linearSolver="@linearSolver")
    solid.addObject('HexahedronSetTopologyContainer', name="topology",
                    src="@../GridTopology/grid")
    solid.addObject('HexahedronSetTopologyModifier')
    dofs = solid.addObject('MechanicalObject', name="dofs",
                           template="Vec3d", src="@topology")
    solid.addObject('LinearSmallStrainFEMForceField', name="FEM",
                    template="Vec3d", youngModulus=E, poissonRatio=nu,
                    topology="@topology")

    eps = 1e-8
    if mode == 'sinus_dirichlet':
        # u_ex=0  
        for name_bc, box in [
            ("xm", [-eps,   -eps,   -eps,    eps,   L+eps, L+eps]),
            ("xp", [L-eps,  -eps,   -eps,   L+eps,  L+eps, L+eps]),
            ("ym", [-eps,   -eps,   -eps,   L+eps,   eps,  L+eps]),
            ("yp", [-eps,  L-eps,   -eps,   L+eps,  L+eps, L+eps]),
            ("zm", [-eps,   -eps,   -eps,   L+eps,  L+eps,  eps  ]),
            ("zp", [-eps,   -eps,  L-eps,   L+eps,  L+eps, L+eps]),
        ]:
            solid.addObject('BoxROI',               name=f"roi_{name_bc}", box=box)
            solid.addObject('FixedProjectiveConstraint',
                            name=f"bc_{name_bc}",
                            indices=f"@roi_{name_bc}.indices")
    else: 
        solid.addObject('BoxROI', name="fixA",
                        box=[-eps,-eps,-eps, eps,eps,eps])
        solid.addObject('FixedProjectiveConstraint', name="bcA",
                        indices="@fixA.indices")
        solid.addObject('BoxROI', name="fixB",
                        box=[L-eps,-eps,-eps, L+eps,eps,eps])
        solid.addObject('PartialFixedProjectiveConstraint', name="bcB",
                        indices="@fixB.indices", fixedDirections=[0,1,1])
        solid.addObject('BoxROI', name="fixC",
                        box=[-eps,L-eps,-eps, eps,L+eps,eps])
        solid.addObject('PartialFixedProjectiveConstraint', name="bcC",
                        indices="@fixC.indices", fixedDirections=[0,0,1])

    # ========================  Forces MMS ==========================================
    n_nodes = nx * ny * nz
    cff = solid.addObject('ConstantForceField', name="MMS_forces",
                          template="Vec3d",
                          indices=list(range(n_nodes)),
                          forces=[[0.,0.,0.]] * n_nodes)

    surf      = solid.addChild('Surface')
    quad_topo = surf.addObject('QuadSetTopologyContainer', name="surfTopo")
    surf.addObject('QuadSetTopologyModifier')
    surf.addObject('Hexa2QuadTopologicalMapping',
                   input="@../topology", output="@surfTopo")
    if visual:
        surf.addObject('OglModel', name="ogl", src="@surfTopo",
                       color=[0.2,0.6,1.0,0.9])
        surf.addObject('IdentityMapping')

    solid.addObject(MMSForceController(
        dofs, cff, quad_topo, nx, ny, nz, E, nu, L, mode,
        name="MMSForceController"))

    return dofs
 

def run_simulation_3d(L, E, nu, nx, ny, nz, mode='linear_neumann'):
    root = Sofa.Core.Node("root")
    root.dt.value = 1.0
    dofs = build_scene_3d(root, L=L, E=E, nu=nu,
                          nx=nx, ny=ny, nz=nz, mode=mode, visual=False)
    Sofa.Simulation.init(root)
    pos_init  = dofs.position.array().copy()
    Sofa.Simulation.animate(root, root.dt.value)
    pos_final = dofs.position.array().copy()
    Sofa.Simulation.unload(root)
    ux = pos_final[:,0] - pos_init[:,0]
    uy = pos_final[:,1] - pos_init[:,1]
    uz = pos_final[:,2] - pos_init[:,2]
    return pos_init, ux, uy, uz

 
 
def l2_error_3d(nodes, ux, uy, uz, L, mode):
    u_ex = get_u_ex(mode)
    ex_x, ex_y, ex_z = u_ex(nodes[:,0], nodes[:,1], nodes[:,2], L)
    return np.sqrt(np.mean((ux-ex_x)**2 + (uy-ex_y)**2 + (uz-ex_z)**2))


def h1_error_3d(nodes, ux, uy, uz, L, nx, ny, nz, mode):
    x, y, z = nodes[:,0], nodes[:,1], nodes[:,2]
    u_ex_fn = get_u_ex(mode)
    ex_x, ex_y, ex_z = u_ex_fn(x, y, z, L)

    h = L / (nx - 1)
    err_u2 = np.sum((ux-ex_x)**2 + (uy-ex_y)**2 + (uz-ex_z)**2) * h**3

    d = get_derivatives(mode)
    all_idx = _build_element_connectivity(nx, ny, nz)  # (n_e, 8)
    xe_all  = nodes[all_idx]                            # (n_e, 8, 3)
    ue_all  = np.column_stack([ux, uy, uz])[all_idx]   # (n_e, 8, 3)

    err_grad2 = 0.0
    for g in range(8):
        J_all              = _dN_GP[g] @ xe_all                   # (n_e, 3, 3)
        detJ_all, Jinv_all = _batch_detinv3(J_all)                # (n_e,), (n_e, 3, 3)

        # dN_phys_all[e,i,j] = ∂N_j/∂x_i
        dN_phys_all = Jinv_all @ _dN_GP[g]                        # (n_e, 3, 8)
        # grad_uh_all[e,i,j] = ∂u_j/∂x_i
        grad_uh_all = dN_phys_all @ ue_all                        # (n_e, 3, 3)

        xg_all = xe_all.transpose(0, 2, 1) @ _N_GP[g]            # (n_e, 3)
        xg0, xg1, xg2 = xg_all[:,0], xg_all[:,1], xg_all[:,2]

        # build (3,3,n_e) then transpose to (n_e,3,3)
        grad_uex_all = np.array([
            [d['dux_dx'](xg0,xg1,xg2,L), d['duy_dx'](xg0,xg1,xg2,L), d['duz_dx'](xg0,xg1,xg2,L)],
            [d['dux_dy'](xg0,xg1,xg2,L), d['duy_dy'](xg0,xg1,xg2,L), d['duz_dy'](xg0,xg1,xg2,L)],
            [d['dux_dz'](xg0,xg1,xg2,L), d['duy_dz'](xg0,xg1,xg2,L), d['duz_dz'](xg0,xg1,xg2,L)],
        ]).transpose(2, 0, 1)  # (n_e, 3, 3)

        err_e = np.sum((grad_uh_all - grad_uex_all)**2, axis=(1, 2))  # (n_e,)
        err_grad2 += np.dot(err_e, np.abs(detJ_all))                  # gw=1

    return np.sqrt(err_u2 + err_grad2)

def plot_convergence(hs, l2s, h1s, mode, out_dir=None):
    if out_dir is None:
        out_dir = RESULTS_DIR
    os.makedirs(out_dir, exist_ok=True)

    hs  = np.array(hs);  l2s = np.array(l2s);  h1s = np.array(h1s)
    ref2 = l2s[-1] * (hs / hs[-1])**2
    ref1 = h1s[-1] * (hs / hs[-1])**1

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    axes[0].loglog(hs, l2s, 'bo-', label='L2 SOFA',   lw=2, markersize=7)
    axes[0].loglog(hs, ref2, 'r--', label='ref $h^2$', lw=1.5)
    axes[0].set(xlabel='h', ylabel='erreur L2',
                title=f'Convergence L2  [{mode}]')
    axes[0].legend(); axes[0].grid(True, which='both', alpha=0.3)
    if len(hs) > 1:
        for k, o in enumerate(np.diff(np.log(l2s))/np.diff(np.log(hs))):
            axes[0].annotate(f'{o:.2f}',
                xy=(np.sqrt(hs[k]*hs[k+1]), np.sqrt(l2s[k]*l2s[k+1])),
                fontsize=9, color='blue')

    axes[1].loglog(hs, h1s, 'gs-', label='H1 SOFA',   lw=2, markersize=7)
    axes[1].loglog(hs, ref1, 'r--', label='ref $h^1$', lw=1.5)
    axes[1].set(xlabel='h', ylabel='erreur H1',
                title=f'Convergence H1  [{mode}]')
    axes[1].legend(); axes[1].grid(True, which='both', alpha=0.3)
    if len(hs) > 1:
        for k, o in enumerate(np.diff(np.log(h1s))/np.diff(np.log(hs))):
            axes[1].annotate(f'{o:.2f}',
                xy=(np.sqrt(hs[k]*hs[k+1]), np.sqrt(h1s[k]*h1s[k+1])),
                fontsize=9, color='green')

    plt.suptitle(f'MMS 3D Hexa Q1 — {mode}', fontsize=13)
    plt.tight_layout()
    out = os.path.join(out_dir, f'convergence_{mode}.png')
    plt.savefig(out, dpi=150); plt.close()
    
 

def plot_displacement(nodes, ux, uy, uz, L, nx, ny, nz, mode, out_dir=None):
    if out_dir is None:
        out_dir = RESULTS_DIR
    os.makedirs(out_dir, exist_ok=True)

    def nidx(i,j,k): return node_idx(i,j,k,nx,ny)
    mid_i=nx//2; mid_j=ny//2; mid_k=nz//2
    x_fine = np.linspace(0, L, 200)
    u_ex_fn = get_u_ex(mode)

    fig, axes = plt.subplots(1, 3, figsize=(15,5))

    line_x = [nidx(i,mid_j,mid_k) for i in range(nx)]
    y_mid=nodes[line_x[0],1]; z_mid=nodes[line_x[0],2]
    ux_ex,_,_ = u_ex_fn(x_fine, np.full_like(x_fine,y_mid), np.full_like(x_fine,z_mid), L)
    axes[0].plot(nodes[line_x,0], ux[line_x], 'bo-', label='SOFA', markersize=7)
    axes[0].plot(x_fine, ux_ex, 'r--', label='MMS exact', lw=2)
    axes[0].set(xlabel='x', ylabel='$u_x$',
                title=f'$u_x$(x)  y={y_mid:.2f}, z={z_mid:.2f}')
    axes[0].legend(); axes[0].grid(alpha=0.3)

    line_y = [nidx(mid_i,j,mid_k) for j in range(ny)]
    x_mid2=nodes[line_y[0],0]; z_mid2=nodes[line_y[0],2]
    _,uy_ex,_ = u_ex_fn(np.full_like(x_fine,x_mid2), x_fine, np.full_like(x_fine,z_mid2), L)
    axes[1].plot(nodes[line_y,1], uy[line_y], 'go-', label='SOFA', markersize=7)
    axes[1].plot(x_fine, uy_ex, 'r--', label='MMS exact', lw=2)
    axes[1].set(xlabel='y', ylabel='$u_y$',
                title=f'$u_y$(y)  x={x_mid2:.2f}, z={z_mid2:.2f}')
    axes[1].legend(); axes[1].grid(alpha=0.3)

    line_z = [nidx(mid_i,mid_j,k) for k in range(nz)]
    x_mid3=nodes[line_z[0],0]; y_mid3=nodes[line_z[0],1]
    _,_,uz_ex = u_ex_fn(np.full_like(x_fine,x_mid3), np.full_like(x_fine,y_mid3), x_fine, L)
    axes[2].plot(nodes[line_z,2], uz[line_z], 'ms-', label='SOFA', markersize=7)
    axes[2].plot(x_fine, uz_ex, 'r--', label='MMS exact', lw=2)
    axes[2].set(xlabel='z', ylabel='$u_z$',
                title=f'$u_z$(z)  x={x_mid3:.2f}, y={y_mid3:.2f}')
    axes[2].legend(); axes[2].grid(alpha=0.3)

    plt.suptitle(f'MMS 3D Hexa Q1 [{mode}] — 1D  ({nx}×{ny}×{nz})', fontsize=13)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir,f'1d_{mode}_nx{nx}.png'), dpi=150); plt.close()

    # Coupes 2D
    fig2, axes2 = plt.subplots(2, 3, figsize=(16,10))

    sl_z = np.array([[nidx(i,j,mid_k) for i in range(nx)] for j in range(ny)])
    X2,Y2 = nodes[sl_z,0], nodes[sl_z,1]
    for row,u,title in [(0,ux,'$u_x$ z=mid'),(1,uy,'$u_y$ z=mid')]:
        im=axes2[row,0].pcolormesh(X2,Y2,u[sl_z],cmap='RdBu_r',shading='auto')
        axes2[row,0].set(xlabel='x',ylabel='y',title=title)
        plt.colorbar(im,ax=axes2[row,0])

    sl_y = np.array([[nidx(i,mid_j,k) for i in range(nx)] for k in range(nz)])
    X2,Z2 = nodes[sl_y,0], nodes[sl_y,2]
    for row,u,title in [(0,ux,'$u_x$ y=mid'),(1,uz,'$u_z$ y=mid')]:
        im=axes2[row,1].pcolormesh(X2,Z2,u[sl_y],cmap='RdBu_r',shading='auto')
        axes2[row,1].set(xlabel='x',ylabel='z',title=title)
        plt.colorbar(im,ax=axes2[row,1])

    sl_x = np.array([[nidx(mid_i,j,k) for j in range(ny)] for k in range(nz)])
    Y2,Z2 = nodes[sl_x,1], nodes[sl_x,2]
    for row,u,title in [(0,uy,'$u_y$ x=mid'),(1,uz,'$u_z$ x=mid')]:
        im=axes2[row,2].pcolormesh(Y2,Z2,u[sl_x],cmap='RdBu_r',shading='auto')
        axes2[row,2].set(xlabel='y',ylabel='z',title=title)
        plt.colorbar(im,ax=axes2[row,2])

    plt.suptitle(f'MMS 3D Hexa Q1 [{mode}] — Displacement ({nx}×{ny}×{nz})', fontsize=13)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir,f'2d_{mode}_nx{nx}.png'), dpi=150); plt.close()




 
def run_convergence_study(mode, nx_list=None, L=1.0, E=1e6, nu=0.3, out_dir=None):
     
    if nx_list is None:
        nx_list = [3, 4, 6, 8, 11, 16]
    
    if out_dir is None:
        out_dir = RESULTS_DIR
    
    os.makedirs(out_dir, exist_ok=True)
    
    hs, l2s, h1s = [], [], []
     
    
    for nx in nx_list:
        ny = nz = nx
        h = L / (nx - 1)
        
        print(f"[{mode}] Simulation nx={nx} (h={h:.4f})...")
        
        # Simulation
        nodes, ux, uy, uz = run_simulation_3d(L, E, nu, nx, ny, nz, mode=mode)
        
        # ============ errors ==============
        l2_err = l2_error_3d(nodes, ux, uy, uz, L, mode)
        h1_err = h1_error_3d(nodes, ux, uy, uz, L, nx, ny, nz, mode)
        
        hs.append(h)
        l2s.append(l2_err)
        h1s.append(h1_err)

        if len(hs) == 1:
            print(f"  L2 error: {l2_err:.6e}")
            print(f"  H1 error: {h1_err:.6e}")
        else:
            l2_rate = np.log(l2s[-2]/l2s[-1]) / np.log(hs[-2]/hs[-1])
            h1_rate = np.log(h1s[-2]/h1s[-1]) / np.log(hs[-2]/hs[-1])
            print(f"  L2 error: {l2_err:.6e}  (rate {l2_rate:.2f})")
            print(f"  H1 error: {h1_err:.6e}  (rate {h1_rate:.2f})")

        if nx == nx_list[-1]:
            plot_displacement(nodes, ux, uy, uz, L, nx, ny, nz, mode, out_dir)

    plot_convergence(hs, l2s, h1s, mode, out_dir)
    
    
    hs_arr = np.array(hs)
    l2_arr = np.array(l2s)
    h1_arr = np.array(h1s)
        
    l2_orders = np.diff(np.log(l2_arr)) / np.diff(np.log(hs_arr))
    h1_orders = np.diff(np.log(h1_arr)) / np.diff(np.log(hs_arr))
        
    print(f"\n  L2 : ordre moyen = {l2_orders.mean():.2f} (attendu: 2.0)")
    print(f"  H1 : ordre moyen = {h1_orders.mean():.2f} (attendu: 1.0)\n")
    
    return {'hs': hs, 'l2s': l2s, 'h1s': h1s}
 

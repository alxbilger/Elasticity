import json
import os
import numpy as np
import matplotlib.pyplot as plt
import Sofa
import Sofa.Core
import Sofa.Simulation

RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results_3d_hexa_linear")

# ======== Parameters mms & derivatives ==================================== 

def load_params(path=None):
    if path is None:
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "params.json")
    with open(path) as f:
        return json.load(f)
 

def lame(E, nu):
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu  = E / (2 * (1 + nu))
    return lam, mu

def u_ex(x, y, z, L): return x / L, y / L, z / L


def stress(E, nu, L):
    lam, mu = lame(E, nu)
    eps = np.diag([1/L, 1/L, 1/L])
    return lam * np.trace(eps) * np.eye(3) + 2 * mu * eps


def traction(x, y, z, nx_c, ny_c, nz_c, E, nu, L):
    t = stress(E, nu, L) @ np.array([nx_c, ny_c, nz_c])
    return t[0], t[1], t[2]
 

def node_idx(i, j, k, nx, ny):  return i + j * nx + k * nx * ny

 

def compute_nodal_forces(nodes, quads, E, nu, L):
    gp = np.array([-1/np.sqrt(3), 1/np.sqrt(3)])
    gw = np.array([1.0, 1.0])
    F  = np.zeros((len(nodes), 3))
    centroid = nodes.mean(axis=0)

    for quad in quads:
        x = nodes[quad]  # (4, 3)
        # determine outward normal sign at quad centre (xi=eta=0)
        dN_dxi0  = np.array([-1.,  1.,  1., -1.]) / 4
        dN_deta0 = np.array([-1., -1.,  1.,  1.]) / 4
        J0   = np.cross(dN_dxi0 @ x, dN_deta0 @ x)
        sign = 1.0 if np.dot(J0, x.mean(axis=0) - centroid) > 0 else -1.0

        for gi, xi in enumerate(gp):
            for gj, eta in enumerate(gp):
                N       = np.array([(1-xi)*(1-eta), (1+xi)*(1-eta),
                                    (1+xi)*(1+eta), (1-xi)*(1+eta)]) / 4
                dN_dxi  = np.array([-(1-eta),  (1-eta),  (1+eta), -(1+eta)]) / 4
                dN_deta = np.array([-(1-xi),  -(1+xi),   (1+xi),  (1-xi)]) / 4
                J_vec   = np.cross(dN_dxi @ x, dN_deta @ x)
                Jmag    = np.linalg.norm(J_vec)
                n_hat   = sign * J_vec / Jmag
                xg      = N @ x
                T       = np.array(traction(*xg, *n_hat, E, nu, L))
                F[quad] += np.outer(N, T) * (gw[gi] * gw[gj] * Jmag)

    res = np.abs(F.sum(axis=0))
    assert res.max() < 1e-6 * np.abs(F).max()
    return F
# =========================== SOFA scene =============================
class MMSForceController(Sofa.Core.Controller):
    def __init__(self, dofs, cff, quad_topo, E, nu, L, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.dofs      = dofs
        self.cff       = cff
        self.quad_topo = quad_topo
        self.E, self.nu, self.L = E, nu, L

    def onSimulationInitDoneEvent(self, event):
        nodes = self.dofs.position.array().copy()
        quads = np.array(self.quad_topo.quads.array())
        F     = compute_nodal_forces(nodes, quads, self.E, self.nu, self.L)
        self.cff.forces.value = F.tolist()


def build_scene_3d(rootNode, L=1.0, E=1e6, nu=0.3, nx=6, ny=6, nz=6,
                   visual=False):
    
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
    ]
    if visual:
        plugins += ["Sofa.Component.Visual", "Sofa.GL.Component.Rendering3D"]

    rootNode.addObject('RequiredPlugin', pluginName=plugins)
    rootNode.addObject('DefaultAnimationLoop')
    if visual:
        rootNode.addObject('DefaultVisualManagerLoop')
        rootNode.addObject('VisualStyle',
                           displayFlags='showVisualModels showWireframe')
    rootNode.gravity.value = [0, 0, 0]
    rootNode.dt.value = 1.0

    # This component has to be added at a different node because one Node
    # cannot have 2 topology components
    regGrid = rootNode.addChild('GridTopology')
    regGrid.addObject('RegularGridTopology',
                    name="grid",
                    nx=nx, ny=ny, nz=nz,
                    min=[0., 0., 0.],
                    max=[L,  L,  L])

    solid = rootNode.addChild('Solid3D')

    solid.addObject('NewtonRaphsonSolver',
                    name="newtonSolver",
                    printLog=False,
                    warnWhenLineSearchFails=True,
                    maxNbIterationsNewton=1,
                    maxNbIterationsLineSearch=1,
                    lineSearchCoefficient=1,
                    relativeSuccessiveStoppingThreshold=0,
                    absoluteResidualStoppingThreshold=1e-7,
                    absoluteEstimateDifferenceThreshold=1e-12,
                    relativeInitialStoppingThreshold=1e-12,
                    relativeEstimateDifferenceThreshold=0)
    solid.addObject('SparseLDLSolver',
                    name="linearSolver",
                    template="CompressedRowSparseMatrixd")
    solid.addObject('StaticSolver',
                    name="staticSolver",
                    newtonSolver="@newtonSolver",
                    linearSolver="@linearSolver")

    solid.addObject('HexahedronSetTopologyContainer',
                    name="topology",
                    src="@../GridTopology/grid")
    solid.addObject('HexahedronSetTopologyModifier')

    dofs = solid.addObject('MechanicalObject',
                           name="dofs",
                           template="Vec3d",
                           src="@topology")

    solid.addObject('LinearSmallStrainFEMForceField',
                    name="FEM",
                    template="Vec3d",
                    youngModulus=E,
                    poissonRatio=nu,
                    topology="@topology")

    #========================== BC using BoxROI ==========================
    eps = 1e-8
    solid.addObject('BoxROI', name="fixA",
                    box=[-eps, -eps, -eps, eps, eps, eps])
    solid.addObject('FixedProjectiveConstraint', name="bcA",
                    indices="@fixA.indices")

    solid.addObject('BoxROI', name="fixB",
                    box=[L-eps, -eps, -eps, L+eps, eps, eps])
    solid.addObject('PartialFixedProjectiveConstraint', name="bcB",
                    indices="@fixB.indices",
                    fixedDirections=[0, 1, 1])

    solid.addObject('BoxROI', name="fixC",
                    box=[-eps, L-eps, -eps, eps, L+eps, eps])
    solid.addObject('PartialFixedProjectiveConstraint', name="bcC",
                    indices="@fixC.indices",
                    fixedDirections=[0, 0, 1])

    n_nodes = nx * ny * nz
    cff = solid.addObject('ConstantForceField',
                          name="MMS_forces",
                          template="Vec3d",
                          indices=list(range(n_nodes)),
                          forces=[[0., 0., 0.]] * n_nodes)

    surf = solid.addChild('Surface')
    quad_topo = surf.addObject('QuadSetTopologyContainer', name="surfTopo")
    surf.addObject('QuadSetTopologyModifier')
    surf.addObject('Hexa2QuadTopologicalMapping',
                   input="@../topology",
                   output="@surfTopo")
    if visual:
        surf.addObject('OglModel',
                       name="ogl",
                       src="@surfTopo",
                       color=[0.2, 0.6, 1.0, 0.9])
        surf.addObject('IdentityMapping')

    solid.addObject(MMSForceController(dofs, cff, quad_topo, E, nu, L,
                                       name="MMSForceController"))

    return dofs

 

def run_simulation_3d(L, E, nu, nx, ny, nz):

    root = Sofa.Core.Node("root")
    root.dt.value = 1.0

    dofs = build_scene_3d(root, L=L, E=E, nu=nu, nx=nx, ny=ny, nz=nz, visual=False)
    Sofa.Simulation.init(root)
    pos_init  = dofs.position.array().copy()
    Sofa.Simulation.animate(root, root.dt.value)
    pos_final = dofs.position.array().copy()
    Sofa.Simulation.unload(root)

    ux = pos_final[:, 0] - pos_init[:, 0]
    uy = pos_final[:, 1] - pos_init[:, 1]
    uz = pos_final[:, 2] - pos_init[:, 2]
    return pos_init, ux, uy, uz

def l2_error_3d(nodes, ux, uy, uz, L):
      
    ex_x, ex_y, ex_z = u_ex(nodes[:, 0], nodes[:, 1], nodes[:, 2], L)
    err = np.sqrt(np.mean((ux - ex_x)**2 + (uy - ex_y)**2 + (uz - ex_z)**2))
    return err

# ======================= Plot ==========================


def plot_displacement(nodes, ux, uy, uz, L, nx, ny, nz):
     
    os.makedirs(RESULTS_DIR, exist_ok=True)
 
    def nidx(i, j, k):
        return node_idx(i, j, k, nx, ny)
 
    mid_i = nx // 2
    mid_j = ny // 2
    mid_k = nz // 2
 
    x_fine = np.linspace(0, L, 200)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
 
    line_x = [nidx(i, mid_j, mid_k) for i in range(nx)]
    xv = nodes[line_x, 0]
    axes[0].plot(xv, ux[line_x], 'bo-', label='SOFA', markersize=7)
    axes[0].plot(x_fine, x_fine / L, 'r--', label='MMS exact', lw=2)
    axes[0].set(xlabel='x', ylabel='$u_x$',
                title=f'$u_x$(x)  y={nodes[line_x[0], 1]:.2f}, z={nodes[line_x[0], 2]:.2f}')
    axes[0].legend(); axes[0].grid(alpha=0.3)
 
    line_y = [nidx(mid_i, j, mid_k) for j in range(ny)]
    yv = nodes[line_y, 1]
    axes[1].plot(yv, uy[line_y], 'go-', label='SOFA', markersize=7)
    axes[1].plot(x_fine, x_fine / L, 'r--', label='MMS exact', lw=2)
    axes[1].set(xlabel='y', ylabel='$u_y$',
                title=f'$u_y$(y)  x={nodes[line_y[0], 0]:.2f}, z={nodes[line_y[0], 2]:.2f}')
    axes[1].legend(); axes[1].grid(alpha=0.3)
 
    line_z = [nidx(mid_i, mid_j, k) for k in range(nz)]
    zv = nodes[line_z, 2]
    axes[2].plot(zv, uz[line_z], 'ms-', label='SOFA', markersize=7)
    axes[2].plot(x_fine, x_fine / L, 'r--', label='MMS exact', lw=2)
    axes[2].set(xlabel='z', ylabel='$u_z$',
                title=f'$u_z$(z)  x={nodes[line_z[0], 0]:.2f}, y={nodes[line_z[0], 1]:.2f}')
    axes[2].legend(); axes[2].grid(alpha=0.3)
 
    plt.suptitle(f'MMS 3D Hexa lin - 1D visual  ({nx}×{ny}×{nz})', fontsize=13)
    plt.tight_layout()
    out1 = os.path.join(RESULTS_DIR, f'1d_visual_nx{nx}.png')
    plt.savefig(out1, dpi=150), plt.close()
    
 
    fig2, axes2 = plt.subplots(2, 3, figsize=(16, 10))
 
    slice_z = np.array([[nidx(i, j, mid_k) for i in range(nx)]
                        for j in range(ny)])         
    X2d  = nodes[slice_z, 0]
    Y2d  = nodes[slice_z, 1]
    UX2d = ux[slice_z]
    UY2d = uy[slice_z]
 
    im = axes2[0, 0].pcolormesh(X2d, Y2d, UX2d, cmap='RdBu_r', shading='auto')
    axes2[0, 0].set(xlabel='x', ylabel='y',
                    title=f'$u_x$  plan z=mid ({nodes[slice_z[0, 0], 2]:.2f})')
    plt.colorbar(im, ax=axes2[0, 0])
 
    im = axes2[1, 0].pcolormesh(X2d, Y2d, UY2d, cmap='RdBu_r', shading='auto')
    axes2[1, 0].set(xlabel='x', ylabel='y',
                    title='$u_y$  plan z=mid')
    plt.colorbar(im, ax=axes2[1, 0])
    slice_y = np.array([[nidx(i, mid_j, k) for i in range(nx)]
                        for k in range(nz)])         
    X2d  = nodes[slice_y, 0]
    Z2d  = nodes[slice_y, 2]
    UX2d = ux[slice_y]
    UZ2d = uz[slice_y]
 
    im = axes2[0, 1].pcolormesh(X2d, Z2d, UX2d, cmap='RdBu_r', shading='auto')
    axes2[0, 1].set(xlabel='x', ylabel='z',
                    title=f'$u_x$  plan y=mid ({nodes[slice_y[0, 0], 1]:.2f})')
    plt.colorbar(im, ax=axes2[0, 1])
 
    im = axes2[1, 1].pcolormesh(X2d, Z2d, UZ2d, cmap='RdBu_r', shading='auto')
    axes2[1, 1].set(xlabel='x', ylabel='z',
                    title='$u_z$  plan y=mid')
    plt.colorbar(im, ax=axes2[1, 1])
 
    slice_x = np.array([[nidx(mid_i, j, k) for j in range(ny)]
                        for k in range(nz)])         
    Y2d  = nodes[slice_x, 1]
    Z2d  = nodes[slice_x, 2]
    UY2d = uy[slice_x]
    UZ2d = uz[slice_x]
 
    im = axes2[0, 2].pcolormesh(Y2d, Z2d, UY2d, cmap='RdBu_r', shading='auto')
    axes2[0, 2].set(xlabel='y', ylabel='z',
                    title=f'$u_y$  plan x=mid ({nodes[slice_x[0, 0], 0]:.2f})')
    plt.colorbar(im, ax=axes2[0, 2])
 
    im = axes2[1, 2].pcolormesh(Y2d, Z2d, UZ2d, cmap='RdBu_r', shading='auto')
    axes2[1, 2].set(xlabel='y', ylabel='z',
                    title='$u_z$  plan x=mid')
    plt.colorbar(im, ax=axes2[1, 2])
 
    plt.suptitle(f'MMS 3D Hexa lin - Displacment  ({nx}×{ny}×{nz})',
                 fontsize=13)
    plt.tight_layout()
    out2 = os.path.join(RESULTS_DIR, f'2d_visual_nx{nx}.png')
    plt.savefig(out2, dpi=150), plt.close()
    

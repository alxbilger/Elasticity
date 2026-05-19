import numpy as np
from utils import (
    load_params,
    run_simulation_3d,
    l2_error_3d,
    h1_error_3d,
    plot_displacement,
    plot_convergence,  
    build_scene_3d,
)

CASE_NAME = "sinus_neumann"
MODE = "sinus_neumann"


def createScene(rootNode):
    cfg = load_params()
    L   = cfg["length"]
    E   = cfg["youngModulus"]
    nu  = cfg["RatioPoisson"]
    nx  = cfg["nx"]
    ny  = cfg["ny"]
    nz  = cfg["nz"]
    build_scene_3d(rootNode, L=L, E=E, nu=nu, nx=nx, ny=ny, nz=nz, 
                   mode=MODE, visual=True)
    return rootNode


if __name__ == "__main__":
    cfg = load_params()
    L   = cfg["length"]
    E   = cfg["youngModulus"]
    nu  = cfg["RatioPoisson"]
    
    # ========== ÉTUDE DE CONVERGENCE ==========
    nx_list = cfg["nxConvergence"][MODE]
    
    hs, l2s, h1s = [], [], []
    
    for nx in nx_list:
        ny = nz = nx
        h = L / (nx - 1)
        
        print(f"\n=== Maillage {nx}×{ny}×{nz} (h={h:.4f}) ===")
        
        nodes, ux, uy, uz = run_simulation_3d(L, E, nu, nx, ny, nz, mode=MODE)
        
        l2_err = l2_error_3d(nodes, ux, uy, uz, L, mode=MODE)
        h1_err = h1_error_3d(nodes, ux, uy, uz, L, nx, ny, nz, mode=MODE)
        
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
    
    
    plot_convergence(hs, l2s, h1s, mode=MODE)
    

    nx = nx_list[-1]
    ny = nz = nx
    nodes, ux, uy, uz = run_simulation_3d(L, E, nu, nx, ny, nz, mode=MODE)
    plot_displacement(nodes, ux, uy, uz, L, nx, ny, nz, mode=MODE)

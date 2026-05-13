import numpy as np

from mms_utils_3d import (
    load_params,
    build_scene_3d,
    run_simulation_3d,
    l2_error_3d,
    plot_displacement,
)

CASE_NAME = "linear3d"

 

def createScene(rootNode):
    cfg = load_params()
    L   = cfg["length"]
    E   = cfg["youngModulus"]
    nu  = cfg["RatioPoisson"]
    nx  = cfg["nx"]
    ny  = cfg["ny"]
    nz  = cfg["nz"]
    build_scene_3d(rootNode, L=L, E=E, nu=nu, nx=nx, ny=ny, nz=nz, visual=True)
    return rootNode
 

if __name__ == "__main__":
    cfg = load_params()
    L   = cfg["length"]
    E   = cfg["youngModulus"]
    nu  = cfg["RatioPoisson"]
    nx  = cfg["nx"]
    ny  = cfg["ny"]
    nz  = cfg["nz"]

    nodes, ux, uy, uz = run_simulation_3d(L, E, nu, nx, ny, nz)

    err = l2_error_3d(nodes, ux, uy, uz, L)
    print(f"L2 nodal error : {err:.4e}")

    plot_displacement(nodes, ux, uy, uz, L, nx, ny, nz)
import json
import os
import sys
import gmsh

RESULTS_DIR = "results_incompressible_q1p1"

def generate_beam2D(length, height, nx, ny, filename):
    """
    Generate a 2D beam mesh (triangular elements)
    """
    gmsh.initialize()
    gmsh.model.add("beam2d")

    # Corner points (tag order matches border traversal in the .edp)
    gmsh.model.geo.addPoint(0,      0,      0, tag=1)
    gmsh.model.geo.addPoint(length, 0,      0, tag=2)
    gmsh.model.geo.addPoint(length, height, 0, tag=3)
    gmsh.model.geo.addPoint(0,      height, 0, tag=4)

    # Boundary lines — tags chosen to match .edp labels directly
    bottom = gmsh.model.geo.addLine(1, 2, tag=2)  # label 2, y=0
    right  = gmsh.model.geo.addLine(2, 3, tag=3)  # label 3, x=length
    top    = gmsh.model.geo.addLine(3, 4, tag=4)  # label 4, y=height
    left   = gmsh.model.geo.addLine(4, 1, tag=1)  # label 1, x=0 (clamped)

    loop = gmsh.model.geo.addCurveLoop([bottom, right, top, left], tag=1)
    surf = gmsh.model.geo.addPlaneSurface([loop], tag=1)

    gmsh.model.geo.synchronize()

    # Physical groups — boundary tags match .edp labels 1–4
    gmsh.model.addPhysicalGroup(1, [left],   tag=1, name="Fixed")
    gmsh.model.addPhysicalGroup(1, [bottom], tag=2, name="Bottom")
    gmsh.model.addPhysicalGroup(1, [right],  tag=3, name="Right")
    gmsh.model.addPhysicalGroup(1, [top],    tag=4, name="Top")
    gmsh.model.addPhysicalGroup(2, [surf],   tag=5, name="Beam")

    # Structured transfinite layout; no recombine → triangular elements
    gmsh.model.mesh.setTransfiniteCurve(bottom, nx)
    gmsh.model.mesh.setTransfiniteCurve(top,    nx)
    gmsh.model.mesh.setTransfiniteCurve(left,   ny)
    gmsh.model.mesh.setTransfiniteCurve(right,  ny)
    gmsh.model.mesh.setTransfiniteSurface(surf)

    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(1)  # ensure linear (first-order) elements only
    _, node_coords, _ = gmsh.model.mesh.getNodes()
    x_positions = node_coords[0::3]
    y_positions = node_coords[1::3]

    os.makedirs(RESULTS_DIR, exist_ok=True)
    msh_path = os.path.join(RESULTS_DIR, filename)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)  # gmshload requires MSH v2
    gmsh.write(msh_path)
    gmsh.finalize()

    return msh_path, x_positions, y_positions


if __name__ == "__main__":
    config_file = sys.argv[1] if len(sys.argv) > 1 else "params.json"
    with open(config_file) as f:
        cfg = json.load(f)

    msh_path, x_positions, y_positions = generate_beam2D(
        length=float(cfg["length"]),
        height=float(cfg["height"]),
        nx=int(cfg["nx"]),
        ny=int(cfg["ny"]),
        filename=cfg.get("meshfile"),
    )
    print("Wrote:", msh_path)
    print("Node x-positions (Gmsh order):", x_positions)
    print("Node y-positions (Gmsh order):", y_positions)

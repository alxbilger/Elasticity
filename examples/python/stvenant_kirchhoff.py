import Sofa
import Elasticity
import Sofa.SofaDeformable
import SofaRuntime
from SofaTypes import Mat3x3

def kronecker(i, j):
    return 1.0 if i == j else 0.0


class Material(Elasticity.HyperElasticMaterialVec3d):
    def __init__(self, young_modulus : float, poisson_ratio : float, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.addData(name="young_modulus", value=young_modulus, type="float", help="Young's modulus", group="Material properties")
        self.addData(name="poisson_ratio", value=poisson_ratio, type="float", help="Poisson's ratio", group="Material properties")

    def secondPiolaKirchhoffStress(self, deformation_gradient):
        F = deformation_gradient
        C = F.transposed() * F
        E = (1.0/2.0) * (C - Mat3x3.Identity())

        lambda_param, mu_param = Sofa.SofaDeformable.toLameParameters3D(self.young_modulus.value,
                                                                        self.poisson_ratio.value)

        return lambda_param * E.trace() * Mat3x3.Identity() + 2.0 * mu_param * E

    def firstPiolaKirchhoffStress(self, deformation_gradient):
        return deformation_gradient * self.secondPiolaKirchhoffStress(deformation_gradient)


    def materialTangentModulus(self, deformation_gradient):

        F = deformation_gradient
        S = self.secondPiolaKirchhoffStress(deformation_gradient)

        lambda_param, mu_param = Sofa.SofaDeformable.toLameParameters3D(self.young_modulus.value,
                                                                        self.poisson_ratio.value)

        def elasticity_tensor(i, j, k, l):
            return mu_param * (kronecker(i, k) * kronecker(j, l) + kronecker(i, l) * kronecker(j, k)) + \
                lambda_param * kronecker(i, j) * kronecker(k, l)

        def tensor_generation(i, j, k, l):
            A_ijkl = kronecker(i, k) * S[l][j]
            for q in range(3):
                for r in range(3):
                    A_ijkl = A_ijkl + F[i][q] * elasticity_tensor(q, j, l, r) * F[k][r]
            return A_ijkl

        return tensor_generation



print(dir(Material))


def createScene(root_node):

    with root_node.addChild("plugins") as plugins:
        plugins.addObject('RequiredPlugin', pluginName="Sofa.Component.Constraint.Projective")
        plugins.addObject('RequiredPlugin', pluginName="Sofa.Component.Engine.Select")
        plugins.addObject('RequiredPlugin', pluginName="Sofa.Component.LinearSolver.Direct")
        plugins.addObject('RequiredPlugin', pluginName="Sofa.Component.LinearSolver.Iterative")
        plugins.addObject('RequiredPlugin', pluginName="Sofa.Component.LinearSolver.Ordering")
        plugins.addObject('RequiredPlugin', pluginName="Sofa.Component.LinearSystem")
        plugins.addObject('RequiredPlugin', pluginName="Sofa.Component.Mass")
        plugins.addObject('RequiredPlugin', pluginName="Sofa.Component.ODESolver.Backward")
        plugins.addObject('RequiredPlugin', pluginName="Sofa.Component.SolidMechanics.FEM.Elastic")
        plugins.addObject('RequiredPlugin', pluginName="Sofa.Component.StateContainer")
        plugins.addObject('RequiredPlugin', pluginName="Sofa.Component.Topology.Container.Dynamic")
        plugins.addObject('RequiredPlugin', pluginName="Sofa.Component.Topology.Container.Grid")
        plugins.addObject('RequiredPlugin', pluginName="Sofa.Component.Topology.Mapping")
        plugins.addObject('RequiredPlugin', pluginName="Sofa.Component.Visual")
        plugins.addObject('RequiredPlugin', pluginName="Sofa.GL.Component.Rendering3D")
        plugins.addObject('RequiredPlugin', pluginName="Elasticity")

    root_node.addObject('DefaultAnimationLoop')

    root_node.addObject('VisualStyle', displayFlags="showBehaviorModels showForceFields")
    root_node.addObject('VisualGrid', size="0.1")
    root_node.addObject('LineAxis', size="0.1")
    root_node.addObject('OglSceneFrame')

    # root_node.addObject('EulerImplicitSolver', name="backward_Euler", rayleighStiffness="0.1", rayleighMass="0.1")
    root_node.addObject('NewtonRaphsonSolver')
    root_node.addObject('StaticSolver', name="solver")
    root_node.addObject('SparseLDLSolver', name="linear_solver", template="CompressedRowSparseMatrixMat3x3d")

    root_node.addObject('RegularGridTopology', name="grid", min="-0.01 -0.01 0", max="0.01 0.01 0.2", n="5 5 30")
    root_node.addObject('MechanicalObject', template="Vec3", name="state", showObject="true")
    root_node.addObject('NodalMassDensity', property="1100")
    root_node.addObject('HexahedronFEMMass')
    root_node.addObject('BoxROI', template="Vec3", name="box_roi", box="-0.011 -0.011 -0.0001   0.011 0.011 0.0001",
                        drawBoxes="1")
    root_node.addObject('FixedProjectiveConstraint', template="Vec3", indices="@box_roi.indices")
    root_node.addObject('HyperelasticityFEMForceField', template="Vec3,Hexahedron", name="fem",
                        computeForceStrategy="sequenced", computeForceDerivStrategy="sequenced")

    young_modulus = 2e5
    poisson_ratio = 0.45
    material = Material(young_modulus, poisson_ratio)
    root_node.addObject(material)

def main():
    import SofaRuntime
    import SofaImGui
    import Sofa.Gui

    root = Sofa.Core.Node("root")
    createScene(root)
    Sofa.Simulation.initRoot(root)

    Sofa.Gui.GUIManager.Init("myscene", "imgui")
    Sofa.Gui.GUIManager.createGUI(root, __file__)
    Sofa.Gui.GUIManager.SetDimension(1600, 900)
    Sofa.Gui.GUIManager.MainLoop(root)
    Sofa.Gui.GUIManager.closeGUI()

if __name__ == "__main__":
    main()
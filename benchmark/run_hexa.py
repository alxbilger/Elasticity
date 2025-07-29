import os
os.environ['BENCHMARK_OUT'] = 'elasticity_hexa.csv'
os.environ['BENCHMARK_OUT_FORMAT'] = 'csv'

import google_benchmark as benchmark
import Sofa
import SofaRuntime
from SofaRuntime import Timer
import Sofa.Gui


# Import all components from the SOFA core
SofaRuntime.importPlugin("Sofa.Component")
# Import the Elasticity plugin
SofaRuntime.importPlugin("Elasticity")


def scene_beam_hexa_assembled_simulation(root, scale, hexahedron_force_field, linear_solver):
    root.addObject('DefaultAnimationLoop')
    root.addObject('DefaultVisualManagerLoop')
    root.addObject('EulerImplicitSolver', name="backward Euler", rayleighStiffness="0.1", rayleighMass="0.1")
    linear_solver(root)
    root.addObject('RegularGridTopology', name="grid", min="-5 -5 0", max="5 5 40", n=[2 * scale, 2 * scale, 10 * scale],)
    root.addObject('MechanicalObject', template="Vec3", name="state", showObject="true")

    root.addObject('MeshMatrixMass', totalMass="1500", topology="@grid")
    hexahedron_force_field(root)

    root.addObject('BoxROI', template="Vec3", name="box_roi", box="-6 -6 -1 6 6 0.1", drawBoxes="1")
    root.addObject('FixedProjectiveConstraint', template="Vec3", indices="@box_roi.indices")

nbIterations = 1

def compute_force_timer(records):
    return records['Simulation::animate']['solve']['Mechanical (root)']['ComputeForce']['total_time']

def compute_addForce_timer(records):
    return records['Simulation::animate']['solve']['Mechanical (root)']['MBKSolve']['CG-Solve']['total_time']

def benchmark_beam(state, hexahedron_force_field, linear_solver, timer):
    while state:
        root = Sofa.Core.Node("root")
        scene_beam_hexa_assembled_simulation(root, state.range(0), hexahedron_force_field, linear_solver)
        Sofa.Simulation.init(root)
        state.counters["nbElements"] = len(root.grid.hexahedra.value)

        if with_gui:
            Sofa.Gui.GUIManager.Init("myscene", "imgui")
            Sofa.Gui.GUIManager.createGUI(root, __file__)
            Sofa.Gui.GUIManager.SetDimension(1080, 1080)
            # Initialization of the scene will be done here
            Sofa.Gui.GUIManager.MainLoop(root)
            Sofa.Gui.GUIManager.closeGUI()

        Timer.setEnabled("Animate", True)
        Timer.setOutputType("Animate", "json")  # hack to avoid printing in the console

        nb_time_steps = 100
        avg_addDForce = 0

        for iteration in range(nb_time_steps):
            Timer.begin("Animate")
            Sofa.Simulation.animate(root, root.dt.value)

            # The first iteration is not added in the measures because precomputations often happen during this iteration
            if iteration != 0:
                records = Timer.getRecords("Animate")
                # print(records)
                addForce_duration = timer(records)
                avg_addDForce += addForce_duration

            Timer.end("Animate")

        Timer.clear()
        Sofa.Simulation.unload(root)

        avg_addDForce = avg_addDForce / (nb_time_steps - 1)
        state.set_iteration_time(avg_addDForce / 1000.)

minScaleFactor = 1
maxScaleFactor = 4

@benchmark.register
@benchmark.option.dense_range(minScaleFactor, maxScaleFactor, 1)
@benchmark.option.unit(benchmark.kMillisecond)
@benchmark.option.use_manual_time()
@benchmark.option.iterations(nbIterations)
def benchmark_beam_hexa_linear_assembled_elasticity_simulation(state):
    def hexahedron_force_field(node):
        node.addObject('LinearSmallStrainFEMForceField', name="FEM", youngModulus="10000",
                        poissonRatio="0.45", topology="@grid", computeVonMisesStress=False)

    def linear_solver(root):
        root.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixd")

    benchmark_beam(state, hexahedron_force_field, linear_solver, compute_force_timer)

@benchmark.register
@benchmark.option.dense_range(minScaleFactor, maxScaleFactor, 1)
@benchmark.option.unit(benchmark.kMillisecond)
@benchmark.option.use_manual_time()
@benchmark.option.iterations(nbIterations)
def benchmark_beam_hexa_linear_assembled_sofa_simulation(state):
    def hexahedron_force_field(node):
        node.addObject('HexahedronFEMForceField', name="FEM", youngModulus="10000", method="small",
                        poissonRatio="0.45", topology="@grid")

    def linear_solver(root):
        root.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixd")

    benchmark_beam(state, hexahedron_force_field, linear_solver, compute_force_timer)

@benchmark.register
@benchmark.option.dense_range(minScaleFactor, maxScaleFactor, 1)
@benchmark.option.unit(benchmark.kMillisecond)
@benchmark.option.use_manual_time()
@benchmark.option.iterations(nbIterations)
def benchmark_beam_hexa_corotational_assembled_elasticity_simulation(state):
    def hexahedron_force_field(node):
        node.addObject('CorotationalFEMForceField', name="FEM", youngModulus="10000",
                        poissonRatio="0.45", topology="@grid", computeVonMisesStress=False)

    def linear_solver(root):
        root.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixd")

    benchmark_beam(state, hexahedron_force_field, linear_solver, compute_force_timer)

@benchmark.register
@benchmark.option.dense_range(minScaleFactor, maxScaleFactor, 1)
@benchmark.option.unit(benchmark.kMillisecond)
@benchmark.option.use_manual_time()
@benchmark.option.iterations(nbIterations)
def benchmark_beam_hexa_corotational_assembled_sofa_simulation(state):
    def hexahedron_force_field(node):
        node.addObject('HexahedronFEMForceField', name="FEM", youngModulus="10000", method="polar",
                        poissonRatio="0.45", topology="@grid")
    def linear_solver(root):
        root.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixd")

    benchmark_beam(state, hexahedron_force_field, linear_solver, compute_force_timer)


@benchmark.register
@benchmark.option.dense_range(minScaleFactor, maxScaleFactor, 1)
@benchmark.option.unit(benchmark.kMillisecond)
@benchmark.option.use_manual_time()
@benchmark.option.iterations(nbIterations)
def benchmark_beam_hexa_corotational_matrixfree_elasticity_simulation(state):
    def hexahedron_force_field(node):
        node.addObject('CorotationalFEMForceField', name="FEM", youngModulus="10000",
                        poissonRatio="0.45", topology="@grid", computeVonMisesStress=False)

    def linear_solver(root):
        root.addObject('CGLinearSolver', iterations=25, tolerance=1.0e-9, threshold=1.0e-9)

    benchmark_beam(state, hexahedron_force_field, linear_solver, compute_addForce_timer)

@benchmark.register
@benchmark.option.dense_range(minScaleFactor, maxScaleFactor, 1)
@benchmark.option.unit(benchmark.kMillisecond)
@benchmark.option.use_manual_time()
@benchmark.option.iterations(nbIterations)
def benchmark_beam_hexa_corotational_matrixfree_sofa_simulation(state):
    def hexahedron_force_field(node):
        node.addObject('HexahedronFEMForceField', name="FEM", youngModulus="10000", method="polar",
                        poissonRatio="0.45", topology="@grid")
    def linear_solver(root):
        root.addObject('CGLinearSolver', iterations=25, tolerance=1.0e-9, threshold=1.0e-9)

    benchmark_beam(state, hexahedron_force_field, linear_solver, compute_addForce_timer)


# Class used to fake a benchmark, so it can run without the Google Benchmark framework
class FakeState:
    def __init__(self, factor):
        self.factor = factor
        self.counters = dict()

    def __bool__(self):
        return True

    def range(self, integer):
        return self.factor

    def set_iteration_time(self, _):
        pass

if __name__ == "__main__":
    with_gui = False

    benchmark.main()

    # The following code is for debugging a SOFA scene using a GUI
    # import SofaImGui
    # with_gui = True
    # benchmark_beam_hexa_corotational_matrixfree_sofa_simulation(FakeState(5))

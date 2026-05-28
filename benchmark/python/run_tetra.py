import os
os.environ['BENCHMARK_OUT'] = 'elasticity_tetra.csv'
os.environ['BENCHMARK_OUT_FORMAT'] = 'csv'

import google_benchmark as benchmark
import Sofa
import SofaRuntime
from SofaRuntime import Timer
import Sofa.Gui


# Import all components from the SOFA core
SofaRuntime.importPlugin("Sofa.Component")


def scene_beam_tetra_assembled_simulation(root, scale, tetrahedron_force_field, linear_solver):
    root.addObject('DefaultAnimationLoop')
    root.addObject('DefaultVisualManagerLoop')
    root.addObject('EulerImplicitSolver', name="backward_euler", rayleighStiffness="0.1", rayleighMass="0.1")
    linear_solver(root)
    root.addObject('RegularGridTopology', name="grid", min="-5 -5 0", max="5 5 40", n=[2 * scale, 2 * scale, 10 * scale],)
    root.addObject('MechanicalObject', template="Vec3", name="state", showObject="true")
    root.addObject('NodalMassDensity', property=0.375)
    tetra = root.addChild('tetra')

    tetra.addObject('TetrahedronSetTopologyContainer', name="Tetra_topo")
    tetra.addObject('TetrahedronSetTopologyModifier', name="Modifier")
    tetra.addObject('TetrahedronSetGeometryAlgorithms', template="Vec3", name="GeomAlgo", drawTetrahedra="true")
    tetra.addObject('Hexa2TetraTopologicalMapping', input="@grid", output="@Tetra_topo")
    tetra.addObject('TetrahedronFEMMass')
    tetrahedron_force_field(tetra)

    root.addObject('BoxROI', template="Vec3", name="box_roi", box="-6 -6 -1 6 6 0.1", drawBoxes="1")
    root.addObject('FixedProjectiveConstraint', template="Vec3", indices="@box_roi.indices")

nbIterations = 1

def compute_force_timer(records):
    return records['Simulation::animate']['solve']['Mechanical (root)']['ComputeForce']['total_time']

def compute_addDForce_timer(records):
    return records['Simulation::animate']['solve']['Mechanical (root)']['MBKSolve']['CG-Solve']['total_time']

def benchmark_beam(state, tetrahedron_force_field, linear_solver, timer):
    while state:
        root = Sofa.Core.Node("root")
        scene_beam_tetra_assembled_simulation(root, state.range(0), tetrahedron_force_field, linear_solver)
        Sofa.Simulation.init(root)
        state.counters["nbElements"] = len(root.grid.hexahedra.value)

        if with_gui:
            Sofa.Gui.GUIManager.Init("myscene", "imgui")
            Sofa.Gui.GUIManager.createGUI(root, __file__)
            Sofa.Gui.GUIManager.SetDimension(1080, 1080)
            # Initialization of the scene will be done here
            Sofa.Gui.GUIManager.MainLoop(root)
            Sofa.Gui.GUIManager.closeGUI()
        else:
            Timer.setEnabled("Animate", True)
            Timer.setOutputType("Animate", "json")  # hack to avoid printing in the console

            nb_time_steps = 10
            avg_measured_duration = 0

            for iteration in range(nb_time_steps):
                Timer.begin("Animate")
                Sofa.Simulation.animate(root, root.dt.value)

                # The first iteration is not added in the measures because precomputations often happen during this iteration
                if iteration != 0:
                    records = Timer.getRecords("Animate")
                    # print(records)
                    measured_duration = timer(records)
                    avg_measured_duration += measured_duration

                Timer.end("Animate")

            Timer.clear()
            Sofa.Simulation.unload(root)

            avg_measured_duration = avg_measured_duration / (nb_time_steps - 1)
            state.set_iteration_time(avg_measured_duration / 1000.)

minScaleFactor = 2
maxScaleFactor = 8
maxScaleFactor_assembled = 5

@benchmark.register
@benchmark.option.dense_range(minScaleFactor, maxScaleFactor_assembled, 1)
@benchmark.option.unit(benchmark.kMillisecond)
@benchmark.option.use_manual_time()
@benchmark.option.iterations(nbIterations)
def benchmark_beam_tetra_linear_assembled_elasticity(state):
    def tetrahedron_force_field(node):
        node.addObject('TetrahedronLinearSmallStrainFEMForceField', name="FEM", youngModulus="10000",
                        poissonRatio="0.45", topology="@Tetra_topo", computeForceStrategy='sequenced', computeForceDerivStrategy='sequenced')

    def linear_solver(root):
        root.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixd")

    benchmark_beam(state, tetrahedron_force_field, linear_solver, compute_force_timer)

@benchmark.register
@benchmark.option.dense_range(minScaleFactor, maxScaleFactor_assembled, 1)
@benchmark.option.unit(benchmark.kMillisecond)
@benchmark.option.use_manual_time()
@benchmark.option.iterations(nbIterations)
def benchmark_beam_tetra_linear_assembled_legacy(state):
    def tetrahedron_force_field(node):
        node.addObject('TetrahedronFEMForceField', name="FEM", youngModulus="10000", method="small",
                        poissonRatio="0.45", topology="@Tetra_topo")

    def linear_solver(root):
        root.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixd")

    benchmark_beam(state, tetrahedron_force_field, linear_solver, compute_force_timer)

@benchmark.register
@benchmark.option.dense_range(minScaleFactor, maxScaleFactor_assembled, 1)
@benchmark.option.unit(benchmark.kMillisecond)
@benchmark.option.use_manual_time()
@benchmark.option.iterations(nbIterations)
def benchmark_beam_tetra_corotational_assembled_elasticity(state):
    def tetrahedron_force_field(node):
        node.addObject('TetrahedronCorotationalFEMForceField', name="FEM", youngModulus="10000",
                        poissonRatio="0.45", topology="@Tetra_topo", computeForceStrategy='sequenced', computeForceDerivStrategy='sequenced')

    def linear_solver(root):
        root.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixd")

    benchmark_beam(state, tetrahedron_force_field, linear_solver, compute_force_timer)

@benchmark.register
@benchmark.option.dense_range(minScaleFactor, maxScaleFactor_assembled, 1)
@benchmark.option.unit(benchmark.kMillisecond)
@benchmark.option.use_manual_time()
@benchmark.option.iterations(nbIterations)
def benchmark_beam_tetra_corotational_assembled_legacy(state):
    def tetrahedron_force_field(node):
        node.addObject('TetrahedronFEMForceField', name="FEM", youngModulus="10000", method="svd",
                        poissonRatio="0.45", topology="@Tetra_topo")
    def linear_solver(root):
        root.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixd")

    benchmark_beam(state, tetrahedron_force_field, linear_solver, compute_force_timer)


@benchmark.register
@benchmark.option.dense_range(minScaleFactor, maxScaleFactor, 1)
@benchmark.option.unit(benchmark.kMillisecond)
@benchmark.option.use_manual_time()
@benchmark.option.iterations(nbIterations)
def benchmark_beam_tetra_corotational_matrixfree_elasticity(state):
    def tetrahedron_force_field(node):
        node.addObject('TetrahedronCorotationalFEMForceField', name="FEM", youngModulus="10000",
                        poissonRatio="0.45", topology="@Tetra_topo", computeForceStrategy='sequenced', computeForceDerivStrategy='sequenced')

    def linear_solver(root):
        root.addObject('CGLinearSolver', iterations=25, tolerance=1.0e-9, threshold=1.0e-9)

    benchmark_beam(state, tetrahedron_force_field, linear_solver, compute_addDForce_timer)

@benchmark.register
@benchmark.option.dense_range(minScaleFactor, maxScaleFactor, 1)
@benchmark.option.unit(benchmark.kMillisecond)
@benchmark.option.use_manual_time()
@benchmark.option.iterations(nbIterations)
def benchmark_beam_tetra_corotational_matrixfree_elasticity_parallel(state):
    def tetrahedron_force_field(node):
        node.addObject('TetrahedronCorotationalFEMForceField', name="FEM", youngModulus="10000",
                       poissonRatio="0.45", topology="@Tetra_topo", computeForceStrategy='parallel', computeForceDerivStrategy='parallel')

    def linear_solver(root):
        root.addObject('CGLinearSolver', iterations=25, tolerance=1.0e-9, threshold=1.0e-9)

    benchmark_beam(state, tetrahedron_force_field, linear_solver, compute_addDForce_timer)

@benchmark.register
@benchmark.option.dense_range(minScaleFactor, maxScaleFactor, 1)
@benchmark.option.unit(benchmark.kMillisecond)
@benchmark.option.use_manual_time()
@benchmark.option.iterations(nbIterations)
def benchmark_beam_tetra_corotational_matrixfree_legacy(state):
    def tetrahedron_force_field(node):
        node.addObject('TetrahedronFEMForceField', name="FEM", youngModulus="10000", method="svd",
                        poissonRatio="0.45", topology="@Tetra_topo")
    def linear_solver(root):
        root.addObject('CGLinearSolver', iterations=25, tolerance=1.0e-9, threshold=1.0e-9)

    benchmark_beam(state, tetrahedron_force_field, linear_solver, compute_addDForce_timer)


# Class used to fake a benchmark, so it can run without the Google Benchmark framework
class FakeState:
    def __init__(self, factor):
        self.factor = factor
        self.counters = dict()
        self.first = True

    def __bool__(self):
        if self.first == True:
            self.first = False
            return True
        else:
            return False

    def range(self, integer):
        return self.factor

    def set_iteration_time(self, _):
        pass

if __name__ == "__main__":
    with_gui = False

    # benchmark.main()

    # The following code is for debugging a SOFA scene using a GUI
    import SofaImGui
    with_gui = True
    benchmark_beam_tetra_corotational_matrixfree_elasticity_parallel(FakeState(3))

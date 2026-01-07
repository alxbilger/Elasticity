#include <Elasticity/component/ElementCorotationalFEMForceField.h>
#include <Elasticity/component/ElementLinearSmallStrainFEMForceField.h>
#include <Elasticity/init.h>
#include <MultiThreading/component/solidmechanics/fem/elastic/ParallelTetrahedronFEMForceField.h>
#include <benchmark/benchmark.h>
#include <sofa/component/init.h>
#include <sofa/component/solidmechanics/fem/elastic/TetrahedronFEMForceField.h>
#include <sofa/component/statecontainer/MechanicalObject.h>
#include <sofa/component/topology/container/dynamic/TetrahedronSetTopologyContainer.h>
#include <sofa/component/topology/container/dynamic/TetrahedronSetTopologyModifier.h>
#include <sofa/component/topology/container/grid/RegularGridTopology.h>
#include <sofa/component/topology/mapping/Hexa2TetraTopologicalMapping.h>
#include <sofa/core/MechanicalParams.h>
#include <sofa/simulation/DefaultAnimationLoop.h>
#include <sofa/simulation/Node.h>
#include <sofa/simulation/Simulation.h>
#include <sofa/simulation/common/init.h>
#include <sofa/simulation/fwd.h>
#include <sofa/simulation/graph/init.h>
#include <sofa/simulation/init.h>

enum class ComponentType : bool
{
    SOFA,
    Elasticity
};

enum class Method : bool
{
    Linear,
    Corotational
};

enum class ComputeStrategy : bool
{
    Parallel,
    Sequential
};

template<ComponentType componentType, Method method, ComputeStrategy computeStrategy>
static sofa::core::behavior::BaseForceField::SPtr createScene(benchmark::State& state)
{
    sofa::simulation::Node::SPtr rootNode = sofa::simulation::getSimulation()->createNewGraph("root");
    rootNode->addObject(sofa::core::objectmodel::New<sofa::simulation::DefaultAnimationLoop>());

    auto hexaTopology = sofa::core::objectmodel::New<sofa::component::topology::container::grid::RegularGridTopology>();
    const auto multiplier = state.range(0);
    hexaTopology->setSize(5 * multiplier, 5 * multiplier, 5 * multiplier);
    hexaTopology->setPos({{0,0,0}, {1,1,1}});
    rootNode->addObject(hexaTopology);

    auto tetraNode = rootNode->createChild("tetra");

    auto tetraTopology = sofa::core::objectmodel::New<sofa::component::topology::container::dynamic::TetrahedronSetTopologyContainer>();
    tetraNode->addObject(tetraTopology);

    auto tetraModifier = sofa::core::objectmodel::New<sofa::component::topology::container::dynamic::TetrahedronSetTopologyModifier>();
    tetraNode->addObject(tetraModifier);

    auto topologicalMapping = sofa::core::objectmodel::New<sofa::component::topology::mapping::Hexa2TetraTopologicalMapping>();
    topologicalMapping->fromModel = hexaTopology;
    topologicalMapping->toModel = tetraTopology;
    tetraNode->addObject(topologicalMapping);

    auto mstate = sofa::core::objectmodel::New<sofa::component::statecontainer::MechanicalObject<sofa::defaulttype::Vec3Types>>();
    auto* topologyPositionData = hexaTopology->findData("position");
    mstate->writeOnlyDx().resize(hexaTopology->getNbPoints());
    mstate->findData("position")->setParent(topologyPositionData);
    tetraNode->addObject(mstate);

    sofa::core::behavior::BaseForceField::SPtr forceField;
    if (componentType == ComponentType::SOFA)
    {
        sofa::component::solidmechanics::fem::elastic::TetrahedronFEMForceField<sofa::defaulttype::Vec3Types>::SPtr tetraForceField;
        if (computeStrategy == ComputeStrategy::Sequential)
        {
            tetraForceField = sofa::core::objectmodel::New<sofa::component::solidmechanics::fem::elastic::TetrahedronFEMForceField<sofa::defaulttype::Vec3Types>>();
        }
        else
        {
            tetraForceField = sofa::core::objectmodel::New<multithreading::component::forcefield::solidmechanics::fem::elastic::ParallelTetrahedronFEMForceField<sofa::defaulttype::Vec3Types>>();
        }

        if (method == Method::Corotational)
        {
            tetraForceField->setMethod("large");
        }
        else
        {
            tetraForceField->setMethod("small");
        }
        tetraForceField->setYoungModulus(1e6);
        tetraForceField->setPoissonRatio(0.45);
        tetraNode->addObject(tetraForceField);
        forceField = tetraForceField;
    }
    else
    {
        if (method == Method::Corotational)
        {
            auto elementForceField = sofa::core::objectmodel::New<elasticity::ElementCorotationalFEMForceField<
               sofa::defaulttype::Vec3Types, sofa::geometry::Tetrahedron>>();
            elementForceField->d_poissonRatio.setValue({0.45});
            elementForceField->d_youngModulus.setValue({1e6});
            {
                auto computeStrategyAccessor = sofa::helper::getWriteOnlyAccessor(elementForceField->d_computeForceStrategy);
                if (computeStrategy == ComputeStrategy::Sequential)
                {
                    computeStrategyAccessor.wref() = elasticity::sequencedComputeStrategy;
                }
                else
                {
                    computeStrategyAccessor.wref() = elasticity::parallelComputeStrategy;
                }
            }
            {
                auto computeStrategyAccessor = sofa::helper::getWriteOnlyAccessor(elementForceField->d_computeForceDerivStrategy);
                if (computeStrategy == ComputeStrategy::Sequential)
                {
                    computeStrategyAccessor.wref() = elasticity::sequencedComputeStrategy;
                }
                else
                {
                    computeStrategyAccessor.wref() = elasticity::parallelComputeStrategy;
                }
            }
            forceField = elementForceField;
        }
        else
        {
            elasticity::ElementLinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Tetrahedron>::SPtr elementForceField
                = sofa::core::objectmodel::New<elasticity::ElementLinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Tetrahedron>>();
            elementForceField->d_poissonRatio.setValue({0.45});
            elementForceField->d_youngModulus.setValue({1e6});
            {
                auto computeStrategyAccessor = sofa::helper::getWriteOnlyAccessor(elementForceField->d_computeForceStrategy);
                if (computeStrategy == ComputeStrategy::Sequential)
                {
                    computeStrategyAccessor.wref() = elasticity::sequencedComputeStrategy;
                }
                else
                {
                    computeStrategyAccessor.wref() = elasticity::parallelComputeStrategy;
                }
            }
            {
                auto computeStrategyAccessor = sofa::helper::getWriteOnlyAccessor(elementForceField->d_computeForceDerivStrategy);
                if (computeStrategy == ComputeStrategy::Sequential)
                {
                    computeStrategyAccessor.wref() = elasticity::sequencedComputeStrategy;
                }
                else
                {
                    computeStrategyAccessor.wref() = elasticity::parallelComputeStrategy;
                }
            }
            forceField = elementForceField;
        }
        tetraNode->addObject(forceField);
    }

    sofa::simulation::node::initRoot(rootNode.get());
    state.counters["NbTetras"] = tetraTopology->getNbTetrahedra();

    return forceField;
}

template<ComponentType componentType, Method method, ComputeStrategy computeStrategy>
static void BM_TetrahedronAddForce(
benchmark::State& state)
{
    auto forcefield = createScene<componentType, method, computeStrategy>(state);

    for (auto _ : state)
    {
        forcefield->addForce(
           sofa::core::MechanicalParams::defaultInstance(), sofa::core::vec_id::write_access::force);
    }
}

template<ComponentType componentType, Method method, ComputeStrategy computeStrategy>
static void BM_TetrahedronAddDForce(
benchmark::State& state)
{
    auto forcefield = createScene<componentType, method, computeStrategy>(state);
    forcefield->addForce( sofa::core::MechanicalParams::defaultInstance(), sofa::core::vec_id::write_access::force);

    for (auto _ : state)
    {
        forcefield->addDForce(
        sofa::core::MechanicalParams::defaultInstance(), sofa::core::vec_id::write_access::dforce);
    }
}

constexpr int minMultiplier = 1;
constexpr int maxMultiplier = 8;
#define BM_OPTIONS ->Unit(benchmark::kMillisecond)->RangeMultiplier(2)->Ranges({ {minMultiplier, maxMultiplier} })

BENCHMARK(BM_TetrahedronAddForce<ComponentType::Elasticity, Method::Linear, ComputeStrategy::Sequential>) BM_OPTIONS;
BENCHMARK(BM_TetrahedronAddForce<ComponentType::SOFA, Method::Linear, ComputeStrategy::Sequential>) BM_OPTIONS;

BENCHMARK(BM_TetrahedronAddForce<ComponentType::Elasticity, Method::Corotational, ComputeStrategy::Sequential>) BM_OPTIONS;
BENCHMARK(BM_TetrahedronAddForce<ComponentType::SOFA, Method::Corotational, ComputeStrategy::Sequential>) BM_OPTIONS;

BENCHMARK(BM_TetrahedronAddForce<ComponentType::Elasticity, Method::Linear, ComputeStrategy::Parallel>) BM_OPTIONS;
BENCHMARK(BM_TetrahedronAddForce<ComponentType::SOFA, Method::Linear, ComputeStrategy::Parallel>) BM_OPTIONS;

BENCHMARK(BM_TetrahedronAddForce<ComponentType::Elasticity, Method::Corotational, ComputeStrategy::Parallel>) BM_OPTIONS;
BENCHMARK(BM_TetrahedronAddForce<ComponentType::SOFA, Method::Corotational, ComputeStrategy::Parallel>) BM_OPTIONS;

BENCHMARK(BM_TetrahedronAddDForce<ComponentType::Elasticity, Method::Linear, ComputeStrategy::Sequential>) BM_OPTIONS;
BENCHMARK(BM_TetrahedronAddDForce<ComponentType::SOFA, Method::Linear, ComputeStrategy::Sequential>) BM_OPTIONS;

BENCHMARK(BM_TetrahedronAddDForce<ComponentType::Elasticity, Method::Corotational, ComputeStrategy::Sequential>) BM_OPTIONS;
BENCHMARK(BM_TetrahedronAddDForce<ComponentType::SOFA, Method::Corotational, ComputeStrategy::Sequential>) BM_OPTIONS;

// BENCHMARK(BM_TetrahedronAddDForce<ComponentType::SOFA, Method::Linear, ComputeStrategy::Parallel>) BM_OPTIONS;
BENCHMARK(BM_TetrahedronAddDForce<ComponentType::Elasticity, Method::Linear, ComputeStrategy::Parallel>) BM_OPTIONS;
// BENCHMARK(BM_TetrahedronAddDForce<ComponentType::SOFA, Method::Corotational, ComputeStrategy::Parallel>) BM_OPTIONS;
BENCHMARK(BM_TetrahedronAddDForce<ComponentType::Elasticity, Method::Corotational, ComputeStrategy::Parallel>) BM_OPTIONS;

int main(int argc, char** argv)
{
    char arg0_default[]="benchmark";
    char* args_default=arg0_default;
    if(!argv)
    {
        argc=1;
        argv=&args_default;
    }
    ::benchmark::Initialize(&argc, argv);
    if(::benchmark::ReportUnrecognizedArguments(argc, argv)) return 1;

    sofa::simulation::graph::init();
    sofa::component::init();
    elasticity::initializePlugin();
    sofa::simulation::core::init();

    ::benchmark::RunSpecifiedBenchmarks();
    ::benchmark::Shutdown();

    sofa::simulation::common::cleanup();
    sofa::simulation::graph::cleanup();
    return 0;
}

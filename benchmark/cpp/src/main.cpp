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

static sofa::core::behavior::BaseForceField::SPtr createScene(
    benchmark::State& state, ComponentType componentType, Method method, ComputeStrategy computeStrategy)
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
            elementForceField->d_poissonRatio.setValue(0.45);
            elementForceField->d_youngModulus.setValue(1e6);
            forceField = elementForceField;
        }
        else
        {
            auto elementForceField = sofa::core::objectmodel::New<elasticity::ElementLinearSmallStrainFEMForceField<
               sofa::defaulttype::Vec3Types, sofa::geometry::Tetrahedron>>();
            elementForceField->d_poissonRatio.setValue(0.45);
            elementForceField->d_youngModulus.setValue(1e6);
            auto computeStrategyAccessor = sofa::helper::getWriteAccessor(elementForceField->d_computeForceStrategy);
            if (computeStrategy == ComputeStrategy::Sequential)
            {
                computeStrategyAccessor->setSelectedItem("sequenced");
            }
            else
            {
                computeStrategyAccessor->setSelectedItem("parallel");
            }
            forceField = elementForceField;
        }
        tetraNode->addObject(forceField);
    }

    sofa::simulation::node::initRoot(rootNode.get());
    state.counters["NbTetras"] = tetraTopology->getNbTetrahedra();

    return forceField;
}

static void BM_TetrahedronAddForce(
benchmark::State& state, ComponentType componentType, Method method, ComputeStrategy computeStrategy)
{
    auto forcefield = createScene(state, componentType, method, computeStrategy);

    for (auto _ : state)
    {
        forcefield->addForce(
           sofa::core::MechanicalParams::defaultInstance(), sofa::core::vec_id::write_access::force);
    }
}

static void BM_TetrahedronAddDForce(
benchmark::State& state, ComponentType componentType, Method method, ComputeStrategy computeStrategy)
{
    auto forcefield = createScene(state, componentType, method, computeStrategy);

    for (auto _ : state)
    {
        forcefield->addDForce(
        sofa::core::MechanicalParams::defaultInstance(), sofa::core::vec_id::write_access::dforce);
    }
}

static void BM_TetrahedronLinearAddForce_SOFASeq(benchmark::State& state)
{
    BM_TetrahedronAddForce(state, ComponentType::SOFA, Method::Linear, ComputeStrategy::Sequential);
}

static void BM_TetrahedronLinearAddDForce_SOFASeq(benchmark::State& state)
{
    BM_TetrahedronAddDForce(state, ComponentType::SOFA, Method::Linear, ComputeStrategy::Sequential);
}



static void BM_TetrahedronLinearAddForce_ElasticitySeq(benchmark::State& state)
{
    BM_TetrahedronAddForce(state, ComponentType::Elasticity, Method::Linear, ComputeStrategy::Sequential);
}

static void BM_TetrahedronLinearAddDForce_ElasticitySeq(benchmark::State& state)
{
    BM_TetrahedronAddDForce(state, ComponentType::Elasticity, Method::Linear, ComputeStrategy::Sequential);
}



static void BM_TetrahedronLinearAddForce_SOFAPar(benchmark::State& state)
{
    BM_TetrahedronAddForce(state, ComponentType::SOFA, Method::Linear, ComputeStrategy::Parallel);
}

static void BM_TetrahedronLinearAddForce_ElasticityPar(benchmark::State& state)
{
    BM_TetrahedronAddForce(state, ComponentType::Elasticity, Method::Linear, ComputeStrategy::Parallel);
}

static void BM_TetrahedronCorotationalAddForce_SOFASeq(benchmark::State& state)
{
    BM_TetrahedronAddForce(state, ComponentType::SOFA, Method::Corotational, ComputeStrategy::Sequential);
}

static void BM_TetrahedronCorotationalAddForce_ElasticitySeq(benchmark::State& state)
{
    BM_TetrahedronAddForce(state, ComponentType::Elasticity, Method::Corotational, ComputeStrategy::Sequential);
}

BENCHMARK(BM_TetrahedronLinearAddForce_SOFASeq)->Unit(benchmark::kMillisecond)->RangeMultiplier(2)->Ranges({ {1, 8} });
BENCHMARK(BM_TetrahedronLinearAddForce_ElasticitySeq)->Unit(benchmark::kMillisecond)->RangeMultiplier(2)->Ranges({ {1, 8} });
BENCHMARK(BM_TetrahedronLinearAddForce_SOFAPar)->Unit(benchmark::kMillisecond)->RangeMultiplier(2)->Ranges({ {1, 8} });
BENCHMARK(BM_TetrahedronLinearAddForce_ElasticityPar)->Unit(benchmark::kMillisecond)->RangeMultiplier(2)->Ranges({ {1, 8} });

BENCHMARK(BM_TetrahedronLinearAddDForce_SOFASeq)->Unit(benchmark::kMillisecond)->RangeMultiplier(2)->Ranges({ {1, 8} });
BENCHMARK(BM_TetrahedronLinearAddDForce_ElasticitySeq)->Unit(benchmark::kMillisecond)->RangeMultiplier(2)->Ranges({ {1, 8} });

BENCHMARK(BM_TetrahedronCorotationalAddForce_SOFASeq)->Unit(benchmark::kMillisecond)->RangeMultiplier(2)->Ranges({ {1, 8} });
BENCHMARK(BM_TetrahedronCorotationalAddForce_ElasticitySeq)->Unit(benchmark::kMillisecond)->RangeMultiplier(2)->Ranges({ {1, 8} });

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

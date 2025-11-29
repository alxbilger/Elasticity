#include <Elasticity/component/ElementCorotationalFEMForceField.h>
#include <Elasticity/init.h>
#include <benchmark/benchmark.h>
#include <sofa/component/init.h>
#include <sofa/component/solidmechanics/fem/elastic/TetrahedronFEMForceField.h>
#include <sofa/component/statecontainer/MechanicalObject.h>
#include <sofa/component/statecontainer/init.h>
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

static void BM_TetrahedronCorotationalAddForce(
    benchmark::State& state, ComponentType componentType, Method method)
{
    sofa::simulation::graph::init();
    sofa::component::init();
    elasticity::initializePlugin();

    sofa::simulation::Node::SPtr rootNode = sofa::simulation::getSimulation()->createNewGraph("root");
    rootNode->addObject(sofa::core::objectmodel::New<sofa::simulation::DefaultAnimationLoop>());

    auto hexaTopology = sofa::core::objectmodel::New<sofa::component::topology::container::grid::RegularGridTopology>();
    hexaTopology->setSize(10, 10, 10);
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
    mstate->findData("position")->setParent(topologyPositionData);
    tetraNode->addObject(mstate);

    sofa::core::behavior::BaseForceField::SPtr forceField;
    if (componentType == ComponentType::SOFA)
    {
        auto tetraForceField = sofa::core::objectmodel::New<sofa::component::solidmechanics::fem::elastic::TetrahedronFEMForceField<sofa::defaulttype::Vec3Types>>();
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
            forceField = elementForceField;
        }
        tetraNode->addObject(forceField);
    }

    sofa::simulation::node::initRoot(rootNode.get());

    for (auto _ : state)
    {
        forceField->addForce(
            sofa::core::MechanicalParams::defaultInstance(), sofa::core::vec_id::write_access::force);
    }

    sofa::simulation::common::cleanup();
    sofa::simulation::graph::cleanup();
}

static void BM_TetrahedronLinearAddForce_SOFA(benchmark::State& state)
{
    BM_TetrahedronCorotationalAddForce(state, ComponentType::SOFA, Method::Linear);
}

static void BM_TetrahedronLinearAddForce_Elasticity(benchmark::State& state)
{
    BM_TetrahedronCorotationalAddForce(state, ComponentType::Elasticity, Method::Linear);
}

static void BM_TetrahedronCorotationalAddForce_SOFA(benchmark::State& state)
{
    BM_TetrahedronCorotationalAddForce(state, ComponentType::SOFA, Method::Corotational);
}

static void BM_TetrahedronCorotationalAddForce_Elasticity(benchmark::State& state)
{
    BM_TetrahedronCorotationalAddForce(state, ComponentType::Elasticity, Method::Corotational);
}

BENCHMARK(BM_TetrahedronLinearAddForce_SOFA)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_TetrahedronLinearAddForce_Elasticity)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_TetrahedronCorotationalAddForce_SOFA)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_TetrahedronCorotationalAddForce_Elasticity)->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();

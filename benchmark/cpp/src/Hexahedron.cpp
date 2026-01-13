
#include <benchmark/benchmark.h>
#include <sofa/component/solidmechanics/fem/elastic/HexahedronFEMForceField.h>
#include <sofa/component/statecontainer/MechanicalObject.h>
#include <sofa/component/topology/container/grid/RegularGridTopology.h>
#include <sofa/core/MechanicalParams.h>
#include <sofa/core/behavior/BaseForceField.h>
#include <sofa/simulation/DefaultAnimationLoop.h>
#include <sofa/simulation/Node.h>
#include <sofa/simulation/Simulation.h>

#include <Elasticity/component/ElementCorotationalFEMForceField.h>
#include <Elasticity/component/ElementLinearSmallStrainFEMForceField.h>

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

template<ComponentType componentType, Method method>
static sofa::core::behavior::BaseForceField::SPtr createScene(benchmark::State& state)
{
    sofa::simulation::Node::SPtr rootNode = sofa::simulation::getSimulation()->createNewGraph("root");
    rootNode->addObject(sofa::core::objectmodel::New<sofa::simulation::DefaultAnimationLoop>());

    auto hexaTopology = sofa::core::objectmodel::New<sofa::component::topology::container::grid::RegularGridTopology>();
    const auto multiplier = state.range(0);
    hexaTopology->setSize(5 * multiplier, 5 * multiplier, 5 * multiplier);
    hexaTopology->setPos({{0,0,0}, {1,1,1}});
    rootNode->addObject(hexaTopology);

    auto mstate = sofa::core::objectmodel::New<sofa::component::statecontainer::MechanicalObject<sofa::defaulttype::Vec3Types>>();
    mstate->writeOnlyDx().resize(hexaTopology->getNbPoints());
    rootNode->addObject(mstate);

    sofa::core::behavior::BaseForceField::SPtr forceField;

    if constexpr (componentType == ComponentType::SOFA)
    {
        auto hexaForceField = sofa::core::objectmodel::New<
            sofa::component::solidmechanics::fem::elastic::HexahedronFEMForceField<
                sofa::defaulttype::Vec3Types>>();
        if constexpr (method == Method::Corotational)
        {
            hexaForceField->d_method.setValue("large");
        }
        else
        {
            hexaForceField->d_method.setValue("small");
        }
        hexaForceField->setYoungModulus(1e6);
        hexaForceField->setPoissonRatio(0.45);
        rootNode->addObject(hexaForceField);
        forceField = hexaForceField;
    }
    else
    {
        if constexpr (method == Method::Corotational)
        {
            auto elementForceField = sofa::core::objectmodel::New<elasticity::ElementCorotationalFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Hexahedron>>();

            elementForceField->d_poissonRatio.setValue({0.45});
            elementForceField->d_youngModulus.setValue({1e6});

            rootNode->addObject(elementForceField);
            forceField = elementForceField;
        }
        else
        {
            auto elementForceField = sofa::core::objectmodel::New<elasticity::ElementLinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Hexahedron>>();

            elementForceField->d_poissonRatio.setValue({0.45});
            elementForceField->d_youngModulus.setValue({1e6});

            rootNode->addObject(elementForceField);
            forceField = elementForceField;
        }
    }

    sofa::simulation::node::initRoot(rootNode.get());
    state.counters["NbHexas"] = hexaTopology->getNbHexahedra();
    return forceField;
}

template<ComponentType componentType, Method method>
static void BM_HexahedronAddForce(
benchmark::State& state)
{
    auto forcefield = createScene<componentType, method>(state);

    for (auto _ : state)
    {
        forcefield->addForce(
           sofa::core::MechanicalParams::defaultInstance(), sofa::core::vec_id::write_access::force);
    }
}

template<ComponentType componentType, Method method>
static void BM_HexahedronAddDForce(
benchmark::State& state)
{
    auto forcefield = createScene<componentType, method>(state);
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

BENCHMARK(BM_HexahedronAddForce<ComponentType::Elasticity, Method::Linear>) BM_OPTIONS;
BENCHMARK(BM_HexahedronAddForce<ComponentType::SOFA, Method::Linear>) BM_OPTIONS;

BENCHMARK(BM_HexahedronAddForce<ComponentType::Elasticity, Method::Corotational>) BM_OPTIONS;
BENCHMARK(BM_HexahedronAddForce<ComponentType::SOFA, Method::Corotational>) BM_OPTIONS;

BENCHMARK(BM_HexahedronAddDForce<ComponentType::Elasticity, Method::Linear>) BM_OPTIONS;
BENCHMARK(BM_HexahedronAddDForce<ComponentType::SOFA, Method::Linear>) BM_OPTIONS;

BENCHMARK(BM_HexahedronAddDForce<ComponentType::Elasticity, Method::Corotational>) BM_OPTIONS;
BENCHMARK(BM_HexahedronAddDForce<ComponentType::SOFA, Method::Corotational>) BM_OPTIONS;

#undef BM_OPTIONS

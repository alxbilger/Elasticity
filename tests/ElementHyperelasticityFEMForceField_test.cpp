#include <Elasticity/component/ElementHyperelasticityFEMForceField.h>
#include <Elasticity/component/material/StVenantKirchhoffMaterial.h>
#include <Elasticity/init.h>
#include <sofa/component/solidmechanics/testing/ForceFieldTestCreation.h>
#include <sofa/component/topology/container/constant/MeshTopology.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/simpleapi/SimpleApi.h>

namespace elasticity
{

template <class DataTypes, class ElementType>
struct ExpectedForce
{
    inline static sofa::VecDeriv_t<DataTypes> value = []()
    {
        sofa::VecDeriv_t<DataTypes> f;
        f.resize(ElementType::NumberOfNodes);
        return f;
    }();
};

template <class DataTypes>
struct ExpectedForce<DataTypes, sofa::geometry::Tetrahedron>
{
    inline static sofa::VecDeriv_t<DataTypes> value = []()
    {
        sofa::VecDeriv_t<DataTypes> f;
        f.resize(sofa::geometry::Tetrahedron::NumberOfNodes);
        constexpr sofa::Real_t<DataTypes> k = 139396.55172413809;
        f[0] = {k, k, k};
        f[1] = {-k, 0, 0};
        f[2] = {0, -k, 0};
        f[3] = {0, 0, -k};
        return f;
    }();
};

template <class DataTypes>
struct ExpectedForce<DataTypes, sofa::geometry::Hexahedron>
{
    inline static sofa::VecDeriv_t<DataTypes> value = []()
    {
        sofa::VecDeriv_t<DataTypes> f;
        f.resize(sofa::geometry::Hexahedron::NumberOfNodes);
        constexpr sofa::Real_t<DataTypes> k = 836379.31034482946;
        f[0] = {k, k, k};
        f[1] = {-k, k, k};
        f[2] = {-k, -k, k};
        f[3] = {k, -k, k};
        f[4] = {k, k, -k};
        f[5] = {-k, k, -k};
        f[6] = {-k, -k, -k};
        f[7] = {k, -k, -k};
        return f;
    }();
};

template <class Real>
struct ExpectedForce<sofa::defaulttype::StdVectorTypes<sofa::type::Vec<3, Real>,sofa::type::Vec<3, Real>, Real>, sofa::geometry::Edge>
{
    using DataTypes = sofa::defaulttype::StdVectorTypes<sofa::type::Vec<3, Real>,sofa::type::Vec<3, Real>, Real>;
    inline static sofa::VecDeriv_t<DataTypes> value = []()
    {
        sofa::VecDeriv_t<DataTypes> f;
        f.resize(sofa::geometry::Edge::NumberOfNodes);
        f[0][0] = 377413.79310344916;
        f[1][0] = -377413.79310344916;
        return f;
    }();
};

/**
 * This test is based on the generic test valid on every force field.
 *
 * It checks the consistency of the derivative functions (addDForce, addKToMatrix,
 * buildStiffnessMatrix) compared to the force function (addForce) using finite differences.
 */
template <class DataTypes, class ElementType>
struct ElementHyperelasticityFEMForceField_stepTest :
    sofa::ForceField_test<ElementHyperelasticityFEMForceField<DataTypes, ElementType>>
{
    using ForceField = ElementHyperelasticityFEMForceField<DataTypes, ElementType>;
    using Inherited = sofa::ForceField_test<ForceField>;
    using DataVecCoord = sofa::DataVecDeriv_t<DataTypes>;
    using DataVecDeriv = sofa::DataVecDeriv_t<DataTypes>;
    using VecCoord = sofa::VecCoord_t<DataTypes>;
    using VecDeriv = sofa::VecDeriv_t<DataTypes>;
    using Coord = sofa::Coord_t<DataTypes>;
    using Deriv = sofa::Deriv_t<DataTypes>;
    using Real = sofa::Real_t<DataTypes>;

    ElementHyperelasticityFEMForceField_stepTest()
    {
        elasticity::initializePlugin();
        sofa::core::ObjectFactory* objectFactory = sofa::core::ObjectFactory::getInstance();
        objectFactory->registerObjectsFromPlugin(elasticity::MODULE_NAME);

        auto topology = sofa::core::objectmodel::New<sofa::component::topology::container::constant::MeshTopology>();
        this->node->addObject(topology);

        //create an element with indices 0, 1, 2, ...
        auto& elementSequence = const_cast<sofa::type::vector<typename FiniteElement<ElementType, DataTypes>::TopologyElement>&>(FiniteElement<ElementType, DataTypes>::getElementSequence(*topology));
        auto& element = elementSequence.emplace_back();
        std::iota(element.begin(), element.end(), 0);

        constexpr auto& referenceElementNodes = FiniteElement<ElementType, DataTypes>::referenceElementNodes;

        auto x = this->dof->writeOnlyRestPositions();
        x.resize(ElementType::NumberOfNodes);
        std::copy(referenceElementNodes.begin(), referenceElementNodes.end(), x.begin());
    }

    void runTestSVK()
    {
        auto material = sofa::core::objectmodel::New<StVenantKirchhoffMaterial<DataTypes>>();
        this->node->addObject(material);

        VecCoord x = this->dof->readRestPositions().ref();
        VecDeriv v,f;

        x.resize(ElementType::NumberOfNodes);
        v.resize(ElementType::NumberOfNodes);
        f.resize(ElementType::NumberOfNodes);

        for (auto& p : x)
        {
            p *= static_cast<Real>(1.1); //extension of 10%
        }

        f = ExpectedForce<DataTypes, ElementType>::value;

        sofa::simulation::node::initRoot(Inherited::node.get());

        Inherited::run_test( x, v, f );
    }
};

template <class TDataTypes, class TElementType>
struct ForceFieldTemplate
{
    using DataTypes = TDataTypes;
    using ElementType = TElementType;
};

/**
 * This alias is introduced because gtest requires only one template parameter
 */
template<class T>
using Hyperelasticity_stepTest =
    ElementHyperelasticityFEMForceField_stepTest<typename T::DataTypes, typename T::ElementType>;

using Tetra3d = ForceFieldTemplate<sofa::defaulttype::Vec3Types, sofa::geometry::Tetrahedron>;
using Hexa3d = ForceFieldTemplate<sofa::defaulttype::Vec3Types, sofa::geometry::Hexahedron>;
using Quad3d = ForceFieldTemplate<sofa::defaulttype::Vec3Types, sofa::geometry::Quad>;
using Quad2d = ForceFieldTemplate<sofa::defaulttype::Vec2Types, sofa::geometry::Quad>;
using Triangle3d = ForceFieldTemplate<sofa::defaulttype::Vec3Types, sofa::geometry::Triangle>;
using Triangle2d = ForceFieldTemplate<sofa::defaulttype::Vec2Types, sofa::geometry::Triangle>;
using Edge3d = ForceFieldTemplate<sofa::defaulttype::Vec3Types, sofa::geometry::Edge>;
using Edge2d = ForceFieldTemplate<sofa::defaulttype::Vec2Types, sofa::geometry::Edge>;
using Edge1d = ForceFieldTemplate<sofa::defaulttype::Vec1Types, sofa::geometry::Edge>;

typedef ::testing::Types<
    Tetra3d,
    Hexa3d,
    Quad3d,
    Quad2d,
    Triangle3d,
    Triangle2d,
    Edge3d,
    Edge2d,
    Edge1d
> TestTypes;

TYPED_TEST_SUITE(Hyperelasticity_stepTest, TestTypes);

TYPED_TEST(Hyperelasticity_stepTest, extension)
{
    this->errorMax *= 1e2;
    this->deltaRange = std::make_pair( 1, this->errorMax * 10 );
    this->debug = true;
    this->flags &= ~sofa::ForceField_test<
        typename Hyperelasticity_stepTest<TypeParam>::ForceField>::TEST_POTENTIAL_ENERGY;

    this->runTestSVK();
}
}

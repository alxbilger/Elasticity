#include <Elasticity/TET4LinearSmallStrainFEMForceField.h>
#include <sofa/component/solidmechanics/testing/ForceFieldTestCreation.h>
#include <sofa/component/topology/container/constant/MeshTopology.h>

namespace sofa::component::topology::container::constant
{
class MeshTopology;
}
namespace elasticity
{

template<class DataTypes>
struct TET4LinearSmallStrainFEMForceField_stepTest : public sofa::ForceField_test<TET4LinearSmallStrainFEMForceField<DataTypes>>
{
    using ForceField = TET4LinearSmallStrainFEMForceField<DataTypes>;
    using Inherited = sofa::ForceField_test<ForceField>;
    using DataVecCoord = sofa::DataVecDeriv_t<DataTypes>;
    using DataVecDeriv = sofa::DataVecDeriv_t<DataTypes>;
    using VecCoord = sofa::VecCoord_t<DataTypes>;
    using VecDeriv = sofa::VecDeriv_t<DataTypes>;
    using Coord = sofa::Coord_t<DataTypes>;
    using Deriv = sofa::Deriv_t<DataTypes>;
    using Real = sofa::Real_t<DataTypes>;
    typedef sofa::type::Vec<3,Real> Vec3;

    TET4LinearSmallStrainFEMForceField_stepTest()
    {
        auto topology = sofa::core::objectmodel::New<sofa::component::topology::container::constant::MeshTopology>();
        this->node->addObject(topology);

        topology->addTetra(0,1,2,3);

        auto x = this->dof->writeOnlyRestPositions();
        x.resize(4);
        DataTypes::set(x[0], 0, 0, 0);
        DataTypes::set(x[1], 1, 0, 0);
        DataTypes::set(x[2], 0, 1, 0);
        DataTypes::set(x[3], 0, 0, 1);
    }

    void runTest()
    {
        VecCoord x = this->dof->readRestPositions().ref();
        VecDeriv v,f;

        //Position
        x.resize(4);
        DataTypes::set(x[0], 0, 0, 0);
        DataTypes::set(x[1], 1.1, 0., 0.);
        DataTypes::set(x[2], 0., 1.1, 0.);
        DataTypes::set(x[3], 0., 0., 1.1);

        //Velocity
        v.resize(4);
        for (auto& vel : v)
        {
            DataTypes::set(vel, 0, 0, 0);
        }

        //mechanical parameters
        this->force->d_poissonRatio.setValue(0);
        this->force->d_youngModulus.setValue(1);

        //Force e*E*S*1/3  = 1*40*sqrt(3)/4*1/3
        f.resize(4);
        constexpr auto k = 1./ 60.;
        DataTypes::set( f[0], k, k, k);
        DataTypes::set( f[1], -k, 0., 0.);
        DataTypes::set( f[2], 0., -k, 0.);
        DataTypes::set( f[3], 0., 0., -k);

        sofa::simulation::node::initRoot(Inherited::node.get());

        Inherited::run_test( x, v, f );
    }
};

typedef ::testing::Types<
    sofa::defaulttype::Vec3Types
> TestTypes;

TYPED_TEST_SUITE(TET4LinearSmallStrainFEMForceField_stepTest, TestTypes);

TYPED_TEST(TET4LinearSmallStrainFEMForceField_stepTest, extension )
{
    this->errorMax *= 1e6;
    this->deltaRange = std::make_pair( 1, this->errorMax * 10 );
    this->debug = true;
    this->flags &= ~sofa::ForceField_test<TET4LinearSmallStrainFEMForceField<TypeParam>>::TEST_POTENTIAL_ENERGY;

    this->runTest();
}

TEST(TET4LinearSmallStrainFEMForceField, computeElasticityTensor)
{
    using Force = TET4LinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types>;

    constexpr auto youngModulus = 1_sreal;
    constexpr auto poissonRatio = 0_sreal;

    const auto C = Force::computeElasticityTensor(youngModulus, poissonRatio);

    for (std::size_t i = 0; i < 3; ++i)
    {
        for (std::size_t j = 0; j < 3; ++j)
        {
            EXPECT_FLOATINGPOINT_EQ(C(i, j), static_cast<SReal>(i == j));
        }
    }

    for (std::size_t i = 3; i < 6; ++i)
    {
        for (std::size_t j = 3; j < 6; ++j)
        {
            EXPECT_FLOATINGPOINT_EQ(C(i, j), static_cast<SReal>(i == j) * 0.5_sreal);
        }
    }

    for (std::size_t i = 3; i < 6; ++i)
    {
        for (std::size_t j = 0; j < 3; ++j)
        {
            EXPECT_FLOATINGPOINT_EQ(C(i, j), 0_sreal);
            EXPECT_FLOATINGPOINT_EQ(C(j, i), 0_sreal);
        }
    }
}

TEST(TET4LinearSmallStrainFEMForceField, computeShapeFunctions)
{
    using Force = TET4LinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types>;

    constexpr std::array<sofa::type::Vec3, 4> tetraNodesCoordinates({
        {0, 0, 0},
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1}
    });

    const auto shapeFunctions = Force::computeShapeFunctions(tetraNodesCoordinates);

    sofa::type::Mat<4, 4, SReal> shapeFunctionMatrix;
    for (sofa::Size i = 0; i < 4; ++i)
    {
        for (sofa::Size j = 0; j < 4; ++j)
        {
            shapeFunctionMatrix(i, j) = shapeFunctions[j][i];
        }
    }

    sofa::type::Mat<4, 4, SReal> coordinatesMatrix;
    for (sofa::Size i = 0; i < 4; ++i)
    {
        coordinatesMatrix(i, 0) = 1;
        for (sofa::Size j = 0; j < 3; ++j)
        {
            coordinatesMatrix(i, j + 1) = tetraNodesCoordinates[i][j];
        }
    }

    const auto expectedIdentity = coordinatesMatrix * shapeFunctionMatrix;
    for (sofa::Size i = 0; i < 4; ++i)
    {
        for (sofa::Size j = 0; j < 4; ++j)
        {
            EXPECT_DOUBLE_EQ(expectedIdentity(i, j), static_cast<SReal>(i == j)) <<  "i = " << i << " j = " << j;
        }
    }
}

TEST(TET4LinearSmallStrainFEMForceField, computeStrainDisplacement)
{
    using Force = TET4LinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types>;

    constexpr std::array<sofa::type::Vec3, 4> tetraNodesCoordinates({
        {0, 0, 0},
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1}
    });

    const auto B = Force::computeStrainDisplacement(tetraNodesCoordinates);

    const sofa::type::Mat<6, 3*4, SReal> expectedB({
        {-1,  0,  0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
        { 0, -1,  0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
        { 0,  0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1},
        {0, -1, -1, 0, 0, 0, 0, 0, 1, 0, 1, 0},
        {-1, 0, -1, 0, 0, 1, 0, 0, 0, 1, 0, 0},
        {-1, -1, 0, 0, 1, 0, 1, 0, 0, 0 ,0, 0}
    });

    for (sofa::Size i = 0; i < 6; ++i)
    {
        for (sofa::Size j = 0; j < 3*4; ++j)
        {
            EXPECT_DOUBLE_EQ(B(i, j), expectedB(i, j)) <<  "i = " << i << " j = " << j;
        }
    }
}

}

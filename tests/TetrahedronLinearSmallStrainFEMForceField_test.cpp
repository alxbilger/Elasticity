#include <Elasticity/component/ElementLinearSmallStrainFEMForceField.h>
#include <Elasticity/finiteelement/FiniteElement[Tetrahedron].h>
#include <Elasticity/impl/LameParameters.h>
#include <Elasticity/impl/MatrixTools.h>
#include <sofa/component/solidmechanics/fem/elastic/TetrahedronFEMForceField.h>
#include <sofa/component/solidmechanics/testing/ForceFieldTestCreation.h>
#include <sofa/component/topology/container/constant/MeshTopology.h>
#include <sofa/core/behavior/BaseForceField.h>
#include <sofa/testing/LinearCongruentialRandomGenerator.h>

namespace elasticity
{

template<class DataTypes>
using TetrahedronLinearSmallStrainFEMForceField =
    ElementLinearSmallStrainFEMForceField<DataTypes, sofa::geometry::Tetrahedron>;

/**
 * This test is based on the generic test valid on every force field.
 *
 * It checks the consistency of the derivative functions (addDForce, addKToMatrix,
 * buildStiffnessMatrix) compared to the force function (addForce) using finite differences.
 */
template <class DataTypes>
struct TET4LinearSmallStrainFEMForceField_stepTest :
    public sofa::ForceField_test<TetrahedronLinearSmallStrainFEMForceField<DataTypes>>
{
    using ForceField = TetrahedronLinearSmallStrainFEMForceField<DataTypes>;
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

        //Position: extension of 10% in each dimension
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
        this->force->d_poissonRatio.setValue({0});
        this->force->d_youngModulus.setValue({1});

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
    this->errorMax *= 1e2;
    this->deltaRange = std::make_pair( 1, this->errorMax * 10 );
    this->debug = true;
    this->flags &= ~sofa::ForceField_test<TetrahedronLinearSmallStrainFEMForceField<TypeParam>>::TEST_POTENTIAL_ENERGY;

    this->runTest();
}

TEST(TET4LinearSmallStrainFEMForceField, computeElasticityTensor)
{
    constexpr auto youngModulus = 1_sreal;
    constexpr auto poissonRatio = 0_sreal;

    const auto [mu, lambda] = toLameParameters<sofa::defaulttype::Vec3Types>(youngModulus, poissonRatio);
    const auto C = makeIsotropicElasticityTensor<sofa::defaulttype::Vec3Types>(mu, lambda).toVoigtMatSym();

    for (std::size_t i = 0; i < 3; ++i)
    {
        for (std::size_t j = 0; j < 3; ++j)
        {
            EXPECT_DOUBLE_EQ(C(i, j), static_cast<SReal>(i == j)) << "i = " << i << " j = " << j;
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

TEST(FiniteElement_Tetra, quadraturePoints)
{
    using FE = FiniteElement<sofa::geometry::Tetrahedron, sofa::defaulttype::Vec3Types>;
    const auto q = FE::quadraturePoints();
    EXPECT_EQ(q.size(), 1);
}

TEST(TET4LinearSmallStrainFEMForceField, jacobian)
{
    using FE = FiniteElement<sofa::geometry::Tetrahedron, sofa::defaulttype::Vec3Types>;

    constexpr std::array<sofa::type::Vec3, 4> tetraNodesCoordinates{{
         {0_sreal, 0_sreal, 0_sreal},
         {1_sreal, 0_sreal, 0_sreal},
         {0_sreal, 1_sreal, 0_sreal},
         {0_sreal, 0_sreal, 1_sreal}
     }};

    const auto q = FE::quadraturePoints();
    const auto dN_dq_ref = FE::gradientShapeFunctions(q[0].first);

    sofa::type::Mat<3, 3, SReal> jacobian;
    for (sofa::Size i = 0; i < 4; ++i)
        jacobian += sofa::type::dyad(tetraNodesCoordinates[i], dN_dq_ref[i]);

    for (sofa::Size i = 0; i < 3; ++i)
        for (sofa::Size j = 0; j < 3; ++j)
           EXPECT_DOUBLE_EQ(jacobian(i, j), static_cast<SReal>(i == j)) << "i = " << i << " j = " << j;

}

struct LegacyComparisonTest : public sofa::testing::BaseSimulationTest, public sofa::testing::NumericTest<SReal>
{
    SceneInstance m_scene{};

    TetrahedronLinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types>::SPtr m_elasticityComponent;
    sofa::component::solidmechanics::fem::elastic::TetrahedronFEMForceField<sofa::defaulttype::Vec3Types>::SPtr m_legacyComponent;
    sofa::component::statecontainer::MechanicalObject<sofa::defaulttype::Vec3Types>::SPtr m_mstate;

    sofa::testing::LinearCongruentialRandomGenerator m_lcrg { 41314325 };

    void doSetUp() override
    {
        m_elasticityComponent = sofa::core::objectmodel::New<TetrahedronLinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types>>();
        m_elasticityComponent->setYoungModulus(1_sreal);
        m_scene.root->addObject(m_elasticityComponent);

        m_legacyComponent = sofa::core::objectmodel::New<sofa::component::solidmechanics::fem::elastic::TetrahedronFEMForceField<sofa::defaulttype::Vec3Types>>();
        m_legacyComponent->d_method.setValue("small");
        m_legacyComponent->setYoungModulus(1_sreal);
        m_scene.root->addObject(m_legacyComponent);

        auto topology = sofa::core::objectmodel::New<sofa::component::topology::container::constant::MeshTopology>();
        topology->addTetra(0,1,2,3);
        m_scene.root->addObject(topology);

        m_mstate = sofa::core::objectmodel::New<sofa::component::statecontainer::MechanicalObject<sofa::defaulttype::Vec3Types>>();
        m_scene.root->addObject(m_mstate);
        m_mstate->resize(4);

        auto x = m_mstate->writeOnlyRestPositions();
        x.resize(4);
        x[0] = { 0, 0, 0 };
        x[1] = { 1, 0, 0 };
        x[2] = { 0, 1, 0 };
        x[3] = { 0, 0, 1 };

        m_mstate->vRealloc(sofa::core::ExecParams::defaultInstance(), sofa::core::vec_id::write_access::dforce);

        m_scene.initScene();
    }

    void resetForces()
    {
        auto forcesAccessor = m_mstate->writeForces();
        forcesAccessor.clear();
    }

    void resetDForces()
    {
        auto dfAccessor = sofa::helper::getWriteAccessor(*m_mstate->write(sofa::core::vec_id::write_access::dforce));
        dfAccessor.clear();
    }

    void checkAddForce(const std::array<sofa::type::Vec3, 4>& vertices, const std::array<sofa::type::Vec3, 4>& vertices_0, const SReal maxError)
    {
        {
            auto x_0 = m_mstate->writeOnlyRestPositions();
            x_0.resize(4);
            for (std::size_t i = 0; i < 4; ++i)
            {
                x_0[i] = vertices_0[i];
            }
        }

        {
            auto x = m_mstate->writeOnlyPositions();
            x.resize(4);
            for (std::size_t i = 0; i < 4; ++i)
            {
                x[i] = vertices[i];
            }
        }

        auto forcesAccessor = m_mstate->writeForces();
        resetForces();

        m_scene.initScene();

        sofa::core::MechanicalParams mparams;
        m_elasticityComponent->toBaseForceField()->addForce(&mparams, sofa::core::vec_id::write_access::force);

        const auto elasticityForces = forcesAccessor.ref();
        resetForces();

        //make sure the forces has been reset
        for (const auto& force : forcesAccessor.ref())
        {
            for (const auto f : force)
            {
                EXPECT_EQ(f, 0_sreal);
            }
        }

        m_legacyComponent->toBaseForceField()->addForce(&mparams, sofa::core::vec_id::write_access::force);
        const auto legacyForces = forcesAccessor.ref();

        for (std::size_t i = 0; i < 4; ++i)
        {
            EXPECT_LT(vectorMaxDiff(legacyForces[i], elasticityForces[i]), maxError)
                << "legacy = " << legacyForces[i] << " vs. elasticity = " << elasticityForces[i];
        }
    }

    void checkAddDForce(const std::array<sofa::type::Vec3, 4>& dx, SReal maxError)
    {
        {
            auto dxAccessor = m_mstate->writeDx();
            dxAccessor.resize(4);
            for (std::size_t i = 0; i < 4; ++i)
            {
                dxAccessor[i] = dx[i];
            }
        }

        sofa::core::MechanicalParams mparams;
        mparams.setKFactor(1._sreal);

        auto dfAccessor = sofa::helper::getWriteAccessor(*m_mstate->write(sofa::core::vec_id::write_access::dforce));

        resetDForces();
        m_elasticityComponent->toBaseForceField()->addDForce(&mparams, sofa::core::vec_id::write_access::dforce);

        const auto elasticityDf = dfAccessor.ref();

        resetDForces();
        m_legacyComponent->toBaseForceField()->addDForce(&mparams, sofa::core::vec_id::write_access::dforce);

        const auto legacyDf = dfAccessor.ref();

        for (std::size_t i = 0; i < 4; ++i)
        {
            EXPECT_LT(vectorMaxDiff(legacyDf[i], elasticityDf[i]), maxError)
                << "legacy = " << legacyDf[i] << " vs. elasticity = " << elasticityDf[i];
        }
    }


    SReal generateScalar(SReal range)
    {
        return m_lcrg.generateInRange(-range, range);
    }

    auto generateVec3(SReal range)
    {
        return sofa::type::Vec3{ generateScalar(range), generateScalar(range), generateScalar(range) };
    }
};

TEST_F(LegacyComparisonTest, checkAddForce)
{
    this->checkAddForce(
        { sofa::type::Vec3{ 0, 0, 0 }, sofa::type::Vec3{ 1, 0, 0 }, sofa::type::Vec3{ 0, 1, 0 }, sofa::type::Vec3{ 0, 0, 1 }},
        { sofa::type::Vec3{ 0, 0, 0 }, sofa::type::Vec3{ 1, 0, 0 }, sofa::type::Vec3{ 0, 1, 0 }, sofa::type::Vec3{ 0, 0, 1 }},
        1e-10_sreal);


    constexpr auto range = 1e4_sreal;

    for (std::size_t i = 0; i < 100; ++i)
    {
        const auto a = generateVec3(range);
        const auto b = generateVec3(range);
        const auto c = generateVec3(range);
        const auto d = generateVec3(range);

        const auto e = generateVec3(range);
        const auto f = generateVec3(range);
        const auto g = generateVec3(range);
        const auto h = generateVec3(range);

        this->checkAddForce({a, b, c, d}, {e, f, g, h}, 1e-3_sreal);
    }
}

TEST_F(LegacyComparisonTest, checkAddDForce)
{
    checkAddDForce({sofa::type::Vec3{1,1,1}, sofa::type::Vec3{1,1,1}, sofa::type::Vec3{1,1,1}, sofa::type::Vec3{1,1,1}}, 1e-8_sreal);

    constexpr auto range = 1e4_sreal;

    for (std::size_t i = 0; i < 100; ++i)
    {
        const auto a = generateVec3(range);
        const auto b = generateVec3(range);
        const auto c = generateVec3(range);
        const auto d = generateVec3(range);

        this->checkAddDForce({a, b, c, d}, 1e-11_sreal);
    }
}

}

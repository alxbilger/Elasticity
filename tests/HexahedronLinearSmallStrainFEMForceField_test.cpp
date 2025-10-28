#include <Elasticity/component/ElementLinearSmallStrainFEMForceField.h>
#include <Elasticity/finiteelement/FiniteElement[Hexahedron].h>
#include <Elasticity/impl/MatrixTools.h>
#include <sofa/component/solidmechanics/testing/ForceFieldTestCreation.h>
#include <sofa/component/topology/container/constant/MeshTopology.h>

namespace elasticity
{
TEST(FiniteElement_Hexa, quadraturePoints)
{
    using FE = FiniteElement<sofa::geometry::Hexahedron, sofa::defaulttype::Vec3Types>;
    const auto q = FE::quadraturePoints();
    EXPECT_EQ(q.size(), 4);
}


TEST(HexahedronLinearSmallStrainFEMForceField, jacobian)
{
    using Force = ElementLinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Hexahedron>;
    using FE = FiniteElement<sofa::geometry::Hexahedron, sofa::defaulttype::Vec3Types>;

    constexpr std::array<sofa::type::Vec3, 8> hexaNodesCoordinates{{
        {-1_sreal, -1_sreal, -1_sreal},
        {1_sreal, -1_sreal, -1_sreal},
        {1_sreal, 1_sreal, -1_sreal},
        {-1_sreal, 1_sreal, -1_sreal},
        {-1_sreal, -1_sreal, 1_sreal},
        {1_sreal, -1_sreal, 1_sreal},
        {1_sreal, 1_sreal, 1_sreal},
        {-1_sreal, 1_sreal, 1_sreal}
     }};

    for (const auto& [q, w] : FE::quadraturePoints())
    {
        const auto dN_dq_ref = FE::gradientShapeFunctions(q);

        sofa::type::Mat<3, 3, SReal> jacobian;
        for (sofa::Size i = 0; i < 8; ++i)
            jacobian += sofa::type::dyad(hexaNodesCoordinates[i], dN_dq_ref[i]);

        for (sofa::Size i = 0; i < 3; ++i)
            for (sofa::Size j = 0; j < 3; ++j)
                EXPECT_NEAR(jacobian(i, j), static_cast<SReal>(i == j), 1e-15) << "i = " << i << " j = " << j << " q = " << q;

        EXPECT_NEAR(jacobian(0, 0), 1., 1e-15);
    }

}

}

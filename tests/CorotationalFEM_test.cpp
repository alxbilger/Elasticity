#include <Elasticity/impl/CorotationalFEM.h>
#include <Elasticity/impl/MatrixTools.h>
#include <Elasticity/finiteelement/FiniteElement[Tetrahedron].h>
#include <Elasticity/finiteelement/FiniteElement[Hexahedron].h>
#include <gtest/gtest.h>

namespace elasticity
{

TEST(CorotationalFEM, centroidTetra)
{
    static const auto centroid = CorotationalFEM<sofa::defaulttype::Vec3Types, sofa::geometry::Tetrahedron>::computeCentroid(
        FiniteElement<sofa::geometry::Tetrahedron, sofa::defaulttype::Vec3Types>::referenceElementNodes);
    EXPECT_EQ(centroid, sofa::type::Vec3(1,1,1)/4);
}
TEST(CorotationalFEM, centroidHexa)
{
    static const auto centroid = CorotationalFEM<sofa::defaulttype::Vec3Types, sofa::geometry::Hexahedron>::computeCentroid(
        FiniteElement<sofa::geometry::Hexahedron, sofa::defaulttype::Vec3Types>::referenceElementNodes);
    EXPECT_EQ(centroid, sofa::type::Vec3());
}

TEST(CorotationalFEM, deformationGradientTetra)
{
    using FE = FiniteElement<sofa::geometry::Tetrahedron, sofa::defaulttype::Vec3Types>;

    constexpr std::array<sofa::type::Vec3, 4> tetraNodesCoordinates{{
         {0_sreal, 0_sreal, 0_sreal},
         {1_sreal, 0_sreal, 0_sreal},
         {0_sreal, 1_sreal, 0_sreal},
         {0_sreal, 0_sreal, 1_sreal}
     }};

    const auto q = FE::quadraturePoints();
    static const auto centroid = CorotationalFEM<sofa::defaulttype::Vec3Types, sofa::geometry::Tetrahedron>::computeCentroid(
        FiniteElement<sofa::geometry::Tetrahedron, sofa::defaulttype::Vec3Types>::referenceElementNodes);
    const auto dN_dq_ref = FE::gradientShapeFunctions(centroid);

    sofa::type::Mat<3, 3, SReal> jacobian;
    for (sofa::Size i = 0; i < 4; ++i)
        jacobian += sofa::type::dyad(tetraNodesCoordinates[i], dN_dq_ref[i]);

    for (sofa::Size i = 0; i < 3; ++i)
        for (sofa::Size j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ(jacobian(i, j), static_cast<SReal>(i == j)) << "i = " << i << " j = " << j;

}

}

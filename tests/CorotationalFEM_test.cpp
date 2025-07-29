#include <Elasticity/impl/CorotationalFEM.h>
#include <Elasticity/impl/MatrixTools.h>
#include <gtest/gtest.h>

namespace elasticity
{

TEST(CorotationalFEM, centroidTetra)
{
    static const auto centroid = CorotationalFEM<sofa::defaulttype::Vec3Types, sofa::geometry::Tetrahedron>::referenceElementCentroid();
    EXPECT_EQ(centroid, sofa::type::Vec3(1,1,1)/4);
}
TEST(CorotationalFEM, centroidHexa)
{
    static const auto centroid = CorotationalFEM<sofa::defaulttype::Vec3Types, sofa::geometry::Hexahedron>::referenceElementCentroid();
    EXPECT_EQ(centroid, sofa::type::Vec3());
}

TEST(CorotationalFEM, deformationGradientTetra)
{
    using FE = FiniteElement<sofa::geometry::Tetrahedron, sofa::defaulttype::Vec3Types>;

    constexpr std::array<sofa::type::Vec3, 4> tetraNodesCoordinates({
         {0, 0, 0},
         {1, 0, 0},
         {0, 1, 0},
         {0, 0, 1}
     });

    const auto q = FE::quadraturePoints();
    const auto dN_dq_ref = FE::gradientShapeFunctions(CorotationalFEM<sofa::defaulttype::Vec3Types, sofa::geometry::Tetrahedron>::referenceElementCentroid());

    const sofa::type::Mat<3, 4, SReal> X_element = nodesMatrix(sofa::topology::Tetrahedron(0,1,2,3), tetraNodesCoordinates);
    const sofa::type::Mat<3, 3, SReal> jacobian = X_element * dN_dq_ref;

    const auto F = CorotationalFEM<sofa::defaulttype::Vec3Types, sofa::geometry::Tetrahedron>::deformationGradient(X_element, jacobian.inverted());

    for (sofa::Size i = 0; i < 3; ++i)
        for (sofa::Size j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ(jacobian(i, j), static_cast<SReal>(i == j)) << "i = " << i << " j = " << j;

}

}

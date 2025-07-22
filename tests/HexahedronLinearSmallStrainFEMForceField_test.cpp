#include <Elasticity/ElementLinearSmallStrainFEMForceField.h>
#include <Elasticity/FiniteElement[Hexahedron].h>
#include <Elasticity/MatrixTools.h>
#include <sofa/component/solidmechanics/testing/ForceFieldTestCreation.h>
#include <sofa/component/topology/container/constant/MeshTopology.h>

namespace elasticity
{
TEST(FiniteElement_Hexa, quadraturePoints)
{
    using FE = FiniteElement<sofa::geometry::Hexahedron, sofa::defaulttype::Vec3Types>;
    const auto q = FE::quadraturePoints();
    EXPECT_EQ(q.size(), 8);
}


TEST(HexahedronLinearSmallStrainFEMForceField, jacobian)
{
    using Force = ElementLinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Hexahedron>;
    using FE = FiniteElement<sofa::geometry::Hexahedron, sofa::defaulttype::Vec3Types>;

    constexpr std::array<sofa::type::Vec3, 8> hexaNodesCoordinates({
        {-1, -1, -1},
        {1, -1, -1},
        {1, 1, -1},
        {-1, 1, -1},
        {-1, -1, 1},
        {1, -1, 1},
        {1, 1, 1},
        {-1, 1, 1}
     });

    const auto X_element = nodesMatrix(sofa::topology::Hexahedron(0,1,2,3,4,5,6,7), hexaNodesCoordinates);
    for (const auto& [q, w] : FE::quadraturePoints())
    {
        const auto dN_dq_ref = FE::gradientShapeFunctions(q);
        const sofa::type::Mat<3, 3, SReal> jacobian = X_element * dN_dq_ref;

        for (sofa::Size i = 0; i < 3; ++i)
            for (sofa::Size j = 0; j < 3; ++j)
                EXPECT_NEAR(jacobian(i, j), static_cast<SReal>(i == j), 1e-15) << "i = " << i << " j = " << j << " q = " << q;

        EXPECT_NEAR(jacobian(0, 0), 1., 1e-15);
    }

}

}

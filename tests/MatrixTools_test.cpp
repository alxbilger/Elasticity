#include <Elasticity/MatrixTools.h>
#include <gtest/gtest.h>
#include <sofa/topology/Tetrahedron.h>
#include <sofa/topology/Triangle.h>

namespace elasticity
{

TEST(Utils, nodesMatrix_triangle)
{
    const auto triangle = sofa::topology::Triangle(0,1,3);
    const sofa::type::vector<sofa::type::Vec3> nodes( {{1,2,3}, {4,5,6}, {7,8,9}, {10,11,12}});
    const auto m = nodesMatrix(triangle, nodes);

    const sofa::type::Mat<3, 3, SReal> expected(
        sofa::type::Vec3({1,4,10}),
        sofa::type::Vec3({2,5,11}),
        sofa::type::Vec3({3,6,12})
        );

    EXPECT_EQ(m, expected);
}

TEST(Utils, nodesMatrix_tetra)
{
    const auto tetra = sofa::topology::Tetrahedron(0,1,3,2);
    const sofa::type::vector<sofa::type::Vec3> nodes( {{1,2,3}, {4,5,6}, {7,8,9}, {10,11,12}});
    const auto m = nodesMatrix(tetra, nodes);

    const sofa::type::Mat<3, 4, SReal> expected(
        sofa::type::Vec4({1,4,10,7}),
        sofa::type::Vec4({2,5,11,8}),
        sofa::type::Vec4({3,6,12,9})
        );

    EXPECT_EQ(m, expected);
}

}

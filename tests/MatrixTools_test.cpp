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

TEST(Utils, pseudoInverse)
{
    constexpr sofa::type::Mat<2, 2, SReal> A = {{-1, 3./2.}, {1., -1.}};
    constexpr sofa::type::Mat<2, 2, SReal> expectedInverse = {{2., 3.}, {2., 2.}};

    const auto pseudoInverse = elasticity::leftPseudoInverse(A);
    const auto inverse = elasticity::inverse(A);

    for (sofa::Size i = 0; i < 2; ++i)
    {
        for (sofa::Size j = 0; j < 2; ++j)
        {
            EXPECT_NEAR(expectedInverse[i][j], inverse[i][j], 1e-15) << "i = " << i << " j = " << j << " A = " << A;
        }
    }

    for (sofa::Size i = 0; i < 2; ++i)
    {
        for (sofa::Size j = 0; j < 2; ++j)
        {
            EXPECT_NEAR(pseudoInverse[i][j], inverse[i][j], 1e-15) << "i = " << i << " j = " << j << " A = " << A;
        }
    }
}

}

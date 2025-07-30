#include <Elasticity/impl/MatrixTools.h>
#include <gtest/gtest.h>
#include <sofa/topology/Tetrahedron.h>
#include <sofa/topology/Triangle.h>
#include <sofa/type/Quat.h>
#include <sofa/testing/LinearCongruentialRandomGenerator.h>

namespace elasticity
{

TEST(MatrixTools, pseudoInverse)
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

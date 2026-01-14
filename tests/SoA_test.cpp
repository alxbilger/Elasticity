#include <Elasticity/impl/MatSoA.h>
#include <gtest/gtest.h>
#include <sofa/testing/LinearCongruentialRandomGenerator.h>

namespace elasticity
{

TEST(SoA, MatrixProduct)
{
    sofa::testing::LinearCongruentialRandomGenerator lcrg(42579813);


    constexpr std::size_t size = 1000;
    elasticity::VecSoA<3, SReal> result(1000);
    elasticity::MatSoA<3, 3, SReal> mat(1000);
    elasticity::VecSoA<3, SReal> vec(1000);

    for (std::size_t i = 0; i < 3; ++i)
    {
        for (auto& e : vec.element(i))
        {
            e = lcrg.generateInUnitRange<SReal>();
        }
    }

    for (std::size_t i = 0; i < 3; ++i)
    {
        for (std::size_t j = 0; j < 3; ++j)
        {
            for (auto& e : mat.element(i, j))
            {
                e = lcrg.generateInUnitRange<SReal>();
            }
        }
    }

    sofa::type::vector<sofa::type::Vec<3, SReal> > vecAoS;
    vec.toAoS(vecAoS);

    sofa::type::vector<sofa::type::Mat<3, 3, SReal> > matAoS;
    mat.toAoS(matAoS);

    elasticity::matrixVectorProduct(result, mat, vec);

    sofa::type::vector<sofa::type::Vec<3, SReal> > resultAoS;
    result.toAoS(resultAoS);

    for (std::size_t s = 0; s < size; ++s)
    {
        const auto expectedResult = matAoS[s] * vecAoS[s];
        for (std::size_t i = 0; i < 3; ++i)
        {
            EXPECT_DOUBLE_EQ(expectedResult[i], resultAoS[s][i]);
        }
    }

}

}

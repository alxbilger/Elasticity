#include <gtest/gtest.h>
#include <sofa/component/solidmechanics/fem/elastic/impl/KroneckerDelta.h>
#include <sofa/core/trait/DataTypes.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/type/FullySymmetric4Tensor.h>

namespace elasticity
{

template<class DataTypes>
void generateThenAccess()
{
    constexpr auto spatial_dimensions = DataTypes::spatial_dimensions;
    const sofa::type::FullySymmetric4Tensor<spatial_dimensions, sofa::Real_t<DataTypes>> t([](sofa::Index i, sofa::Index j, sofa::Index k, sofa::Index l)
    {
        return sofa::component::solidmechanics::fem::elastic::kroneckerDelta<SReal>(i, j) *
                   sofa::component::solidmechanics::fem::elastic::kroneckerDelta<SReal>(k, l);
    });

    for (sofa::Index i = 0; i < spatial_dimensions; ++i)
    {
        for (sofa::Index j = 0; j < spatial_dimensions; ++j)
        {
            for (sofa::Index k = 0; k < spatial_dimensions; ++k)
            {
                for (sofa::Index l = 0; l < spatial_dimensions; ++l)
                {
                    EXPECT_DOUBLE_EQ(t(i, j, k, l), static_cast<SReal>(i == j) * static_cast<SReal>(k == l)) << "i = " << i << " j = " << j << " k = " << k << " l = " << l;
                }
            }
        }
    }
}

TEST(ElasticityTensor, generateThenAccess1d)
{
    generateThenAccess<sofa::defaulttype::Vec1Types>();
}
TEST(ElasticityTensor, generateThenAccess2d)
{
    generateThenAccess<sofa::defaulttype::Vec2Types>();
}
TEST(ElasticityTensor, generateThenAccess3d)
{
    generateThenAccess<sofa::defaulttype::Vec3Types>();
}

}

#pragma once

#include <Elasticity/impl/SymmetricTensor.h>
#include <sofa/core/trait/DataTypes.h>
#include <sofa/type/Mat.h>

namespace elasticity
{

template<class DataTypes>
std::pair<sofa::Real_t<DataTypes>, sofa::Real_t<DataTypes>>
toLameParameters(sofa::Real_t<DataTypes> youngModulus, sofa::Real_t<DataTypes> poissonRatio)
{
    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;
    const auto mu = youngModulus / (2 * (1 + poissonRatio));
    const auto lambda = youngModulus * poissonRatio / ((1 + poissonRatio) * (1 - (spatial_dimensions - 1) * poissonRatio));
    return std::make_pair(mu, lambda);
}

template <class DataTypes>
using ElasticityTensor = sofa::type::Mat<
    symmetric_tensor::NumberOfIndependentElements<DataTypes::spatial_dimensions>,
    symmetric_tensor::NumberOfIndependentElements<DataTypes::spatial_dimensions>,
    sofa::Real_t<DataTypes>>;

template <class DataTypes>
ElasticityTensor<DataTypes> makeIsotropicElasticityTensor(sofa::Real_t<DataTypes> youngModulus, sofa::Real_t<DataTypes> poissonRatio)
{
    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;
    static constexpr sofa::Size NumberOfIndependentElements = symmetric_tensor::NumberOfIndependentElements<spatial_dimensions>;
    using Real = sofa::Real_t<DataTypes>;

    const auto [mu, lambda] = toLameParameters<DataTypes>(youngModulus, poissonRatio);
    ElasticityTensor<DataTypes> C;

    for (sofa::Size i = 0; i < spatial_dimensions; ++i)
    {
        for (sofa::Size j = 0; j < spatial_dimensions; ++j)
        {
            C(i, j) = lambda;
        }
        C(i, i) += 2 * mu;
    }

    for (sofa::Size i = spatial_dimensions; i < NumberOfIndependentElements; ++i)
    {
        C(i, i) = mu;
    }

    return C;
}

}

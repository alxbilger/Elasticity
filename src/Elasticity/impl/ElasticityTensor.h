#pragma once

#include <Elasticity/impl/ComputationStrategy.h>
#include <Elasticity/impl/SymmetricTensor.h>
#include <sofa/core/trait/DataTypes.h>
#include <sofa/type/Mat.h>
#include <Elasticity/impl/StaticSparseMatrix.h>

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
using ElasticityTensorDense = sofa::type::Mat<
    symmetric_tensor::NumberOfIndependentElements<DataTypes::spatial_dimensions>,
    symmetric_tensor::NumberOfIndependentElements<DataTypes::spatial_dimensions>,
    sofa::Real_t<DataTypes>>;


template <class DataTypes>
constexpr auto NumberNonZerosElasticityTensor =
    DataTypes::spatial_dimensions * DataTypes::spatial_dimensions +
    symmetric_tensor::NumberOfIndependentElements<DataTypes::spatial_dimensions> - DataTypes::spatial_dimensions;


template <class DataTypes>
consteval std::array<PairIndices, NumberNonZerosElasticityTensor<DataTypes>> computeListNonZeroEntriesInElasticityTensor()
{
    std::array<PairIndices, NumberNonZerosElasticityTensor<DataTypes>> pairs;
    sofa::Size counter{};
    for (sofa::Size i = 0; i < DataTypes::spatial_dimensions; ++i)
    {
        for (sofa::Size j = 0; j < DataTypes::spatial_dimensions; ++j)
        {
            pairs[counter++] = std::make_pair(i, j);
        }
    }
    for (sofa::Size i = DataTypes::spatial_dimensions; i < symmetric_tensor::NumberOfIndependentElements<DataTypes::spatial_dimensions>; ++i)
    {
        pairs[counter++] = std::make_pair(i, i);
    }
    return pairs;
}

template <class DataTypes>
using ElasticityTensorSparse = StaticSparseMatrix<
    symmetric_tensor::NumberOfIndependentElements<DataTypes::spatial_dimensions>,
    symmetric_tensor::NumberOfIndependentElements<DataTypes::spatial_dimensions>,
    NumberNonZerosElasticityTensor<DataTypes>,
    computeListNonZeroEntriesInElasticityTensor<DataTypes>(),
    sofa::Real_t<DataTypes>
>;

template <class DataTypes, ComputationStrategy strategy>
using ElasticityTensor = std::conditional_t<strategy == ComputationStrategy::SPARSE,
    ElasticityTensorSparse<DataTypes>,
    ElasticityTensorDense<DataTypes>>;

template <class DataTypes, ComputationStrategy strategy>
ElasticityTensor<DataTypes, strategy> makeIsotropicElasticityTensor(sofa::Real_t<DataTypes> youngModulus, sofa::Real_t<DataTypes> poissonRatio)
{
    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;
    static constexpr sofa::Size NumberOfIndependentElements = symmetric_tensor::NumberOfIndependentElements<spatial_dimensions>;
    using Real = sofa::Real_t<DataTypes>;

    const auto [mu, lambda] = toLameParameters<DataTypes>(youngModulus, poissonRatio);
    ElasticityTensor<DataTypes, strategy> C;

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

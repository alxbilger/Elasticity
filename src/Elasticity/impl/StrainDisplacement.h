#pragma once
#include <sofa/core/trait/DataTypes.h>
#include <sofa/type/Mat.h>

#include <Elasticity/impl/SymmetricTensor.h>

namespace elasticity
{

template <class DataTypes, class ElementType>
using StrainDisplacementDense = sofa::type::Mat<
    symmetric_tensor::NumberOfIndependentElements<DataTypes::spatial_dimensions>,
    ElementType::NumberOfNodes * DataTypes::spatial_dimensions,
    sofa::Real_t<DataTypes>>;

template <class DataTypes, class ElementType>
constexpr auto NumberNonZerosStrainDisplacementTensor =
    ElementType::NumberOfNodes * (DataTypes::spatial_dimensions * DataTypes::spatial_dimensions);

template <class DataTypes, class ElementType>
consteval std::array<PairIndices, NumberNonZerosStrainDisplacementTensor<DataTypes, ElementType>> computeListNonZeroEntriesInStrainDisplacementTensor()
{
    std::array<PairIndices, NumberNonZerosStrainDisplacementTensor<DataTypes, ElementType>> pairs;
    sofa::Size counter{};

    for (sofa::Size ne = 0; ne < ElementType::NumberOfNodes; ++ne)
    {
        for (sofa::Size i = 0; i < DataTypes::spatial_dimensions; ++i)
        {
            pairs[counter++] = std::make_pair(i, ne * DataTypes::spatial_dimensions + i);
        }

        auto row = DataTypes::spatial_dimensions;
        for (sofa::Size i = 0; i < DataTypes::spatial_dimensions; ++i)
        {
            for (sofa::Size j = i + 1; j < DataTypes::spatial_dimensions; ++j)
            {
                pairs[counter++] = std::make_pair(row, ne * DataTypes::spatial_dimensions + i);
                pairs[counter++] = std::make_pair(row, ne * DataTypes::spatial_dimensions + j);
                ++row;
            }
        }
    }

    return pairs;
}

template <class DataTypes, class ElementType>
using StrainDisplacementSparse = StaticSparseMatrix<
    symmetric_tensor::NumberOfIndependentElements<DataTypes::spatial_dimensions>,
    ElementType::NumberOfNodes * DataTypes::spatial_dimensions,
    NumberNonZerosStrainDisplacementTensor<DataTypes, ElementType>,
    computeListNonZeroEntriesInStrainDisplacementTensor<DataTypes, ElementType>(),
    sofa::Real_t<DataTypes>>;

template <class DataTypes, class ElementType, ComputationStrategy strategy>
using StrainDisplacement = std::conditional_t<strategy == ComputationStrategy::SPARSE,
    StrainDisplacementSparse<DataTypes, ElementType>,
    StrainDisplacementDense<DataTypes, ElementType>>;

template <class DataTypes, class ElementType, ComputationStrategy strategy>
StrainDisplacement<DataTypes, ElementType, strategy> makeStrainDisplacement(
    const sofa::type::Mat<ElementType::NumberOfNodes, DataTypes::spatial_dimensions, sofa::Real_t<DataTypes> > gradientShapeFunctions)
{
    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;
    static constexpr sofa::Size NumberOfNodesInElement = ElementType::NumberOfNodes;

    StrainDisplacement<DataTypes, ElementType, strategy> B;
    for (sofa::Size ne = 0; ne < NumberOfNodesInElement; ++ne)
    {
        for (sofa::Size i = 0; i < spatial_dimensions; ++i)
        {
            B(i, ne * spatial_dimensions + i) = gradientShapeFunctions[ne][i];
        }

        auto row = spatial_dimensions;
        for (sofa::Size i = 0; i < spatial_dimensions; ++i)
        {
            for (sofa::Size j = i + 1; j < spatial_dimensions; ++j)
            {
                B(row, ne * spatial_dimensions + i) = gradientShapeFunctions[ne][j];
                B(row, ne * spatial_dimensions + j) = gradientShapeFunctions[ne][i];
                ++row;
            }
        }
    }

    return B;
}

}

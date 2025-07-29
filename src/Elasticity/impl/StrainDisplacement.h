#pragma once
#include <sofa/core/trait/DataTypes.h>
#include <sofa/type/Mat.h>

#include <Elasticity/impl/SymmetricTensor.h>

namespace elasticity
{

template <class DataTypes, class ElementType>
using StrainDisplacement = sofa::type::Mat<
    symmetric_tensor::NumberOfIndependentElements<DataTypes::spatial_dimensions>,
    ElementType::NumberOfNodes * DataTypes::spatial_dimensions,
    sofa::Real_t<DataTypes>>;

template <class DataTypes, class ElementType>
StrainDisplacement<DataTypes, ElementType> makeStrainDisplacement(
    const sofa::type::Mat<ElementType::NumberOfNodes, DataTypes::spatial_dimensions, sofa::Real_t<DataTypes> > gradientShapeFunctions)
{
    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;
    static constexpr sofa::Size NumberOfNodesInElement = ElementType::NumberOfNodes;

    StrainDisplacement<DataTypes, ElementType> B;
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

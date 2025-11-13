#pragma once

#include <Elasticity/finiteelement/FiniteElement.h>
#include <Elasticity/impl/FullySymmetric4Tensor.h>
#include <Elasticity/impl/MatrixTools.h>
#include <Elasticity/impl/StrainDisplacement.h>
#include <sofa/type/Mat.h>

namespace elasticity
{

template <class DataTypes, class ElementType>
using ElementStiffness = sofa::type::Mat<
    ElementType::NumberOfNodes * DataTypes::spatial_dimensions,
    ElementType::NumberOfNodes * DataTypes::spatial_dimensions,
    sofa::Real_t<DataTypes>
>;

template <class DataTypes, class ElementType>
ElementStiffness<DataTypes, ElementType> integrate(
    const std::array<sofa::Coord_t<DataTypes>, ElementType::NumberOfNodes>& nodesCoordinates,
    const FullySymmetric4Tensor<DataTypes>& elasticityTensor)
{
    using Real = sofa::Real_t<DataTypes>;
    using FiniteElement = FiniteElement<ElementType, DataTypes>;

    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;
    static constexpr sofa::Size NumberOfNodesInElement = ElementType::NumberOfNodes;
    static constexpr sofa::Size ElementDimension = FiniteElement::ElementDimension;

    ElementStiffness<DataTypes, ElementType> K;

    for (const auto& [quadraturePoint, weight] : FiniteElement::quadraturePoints())
    {
        // gradient of shape functions in the reference element evaluated at the quadrature point
        const sofa::type::Mat<NumberOfNodesInElement, ElementDimension, Real> dN_dq_ref =
            FiniteElement::gradientShapeFunctions(quadraturePoint);

        // jacobian of the mapping from the reference space to the physical space, evaluated at the
        // quadrature point
        sofa::type::Mat<spatial_dimensions, ElementDimension, Real> jacobian;
        for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
            jacobian += sofa::type::dyad(nodesCoordinates[i], dN_dq_ref[i]);

        const auto detJ = elasticity::determinant(jacobian);
        const sofa::type::Mat<ElementDimension, spatial_dimensions, Real> J_inv =
            elasticity::inverse(jacobian);

        // gradient of the shape functions in the physical element evaluated at the quadrature point
        sofa::type::Mat<NumberOfNodesInElement, spatial_dimensions, Real> dN_dq(sofa::type::NOINIT);
        for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
            dN_dq[i] = J_inv.transposed() * dN_dq_ref[i];

        const auto B = makeStrainDisplacement<DataTypes, ElementType>(dN_dq);
        K += (weight * detJ) * B.transposed() * (elasticityTensor.toVoigtMatSym() * B);
    }
    return K;
}



}

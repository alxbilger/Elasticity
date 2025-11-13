#pragma once
#include <sofa/core/topology/BaseMeshTopology.h>

namespace elasticity
{

template <class ElementType, class DataTypes>
struct FiniteElement;

#define FINITEELEMENT_HEADER(ElType, DataTypes, dimension) \
    using Coord = sofa::Coord_t<DataTypes>;\
    using Real = sofa::Real_t<DataTypes>;\
    using ElementType = ElType;\
    using TopologyElement = sofa::topology::Element<ElementType>;\
    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;\
    static constexpr sofa::Size NumberOfNodesInElement = ElementType::NumberOfNodes;\
    static constexpr sofa::Size ElementDimension = dimension;\
    using ReferenceCoord = sofa::type::Vec<ElementDimension, Real>;\
    using ShapeFunctionType = std::function<Real(const ReferenceCoord&)>;\
    using QuadraturePoint = ReferenceCoord; \
    using QuadraturePointAndWeight = std::pair<QuadraturePoint, Real>


template <class ElementType, class DataTypes>
constexpr auto gradientShapeFunctionAtQuadraturePoints()
{
    using FiniteElement = elasticity::FiniteElement<ElementType, DataTypes>;
    using Gradient = sofa::type::Mat<FiniteElement::NumberOfNodesInElement, FiniteElement::ElementDimension, typename FiniteElement::Real>;

    constexpr auto quadraturePoints = FiniteElement::quadraturePoints();

    std::array<Gradient, quadraturePoints.size()> gradients;
    for (sofa::Size i = 0; i < quadraturePoints.size(); ++i)
    {
        gradients[i] = FiniteElement::gradientShapeFunctions(quadraturePoints[i].first);
    }
    return gradients;
};

}

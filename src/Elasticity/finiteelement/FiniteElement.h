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

}

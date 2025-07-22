#pragma once
#include <Elasticity/FiniteElement.h>

namespace elasticity
{

template <class DataTypes>
struct FiniteElement<sofa::geometry::Quad, DataTypes>
{
    FINITEELEMENT_HEADER(sofa::geometry::Quad, DataTypes, 2);
    static_assert(spatial_dimensions > 1, "Quads cannot be defined in 1D");

    constexpr static Real volume(const std::array<Coord, NumberOfNodesInElement>& nodesCoordinates)
    {
        return ElementType::area(nodesCoordinates[0], nodesCoordinates[1], nodesCoordinates[2], nodesCoordinates[3]);
    }

    constexpr static std::array<Coord, NumberOfNodesInElement> referenceElementNodes = []()
    {
        std::array<Coord, NumberOfNodesInElement> nodes;
        DataTypes::set(nodes[0], -1, -1, 0);
        DataTypes::set(nodes[1], 1, -1, 0);
        DataTypes::set(nodes[2], 1, 1, 0);
        DataTypes::set(nodes[3], -1, 1, 0);
        return nodes;
    }();

    static sofa::type::vector<TopologyElement> getElementSequence(sofa::core::topology::BaseMeshTopology& topology)
    {
        return topology.getQuads();
    }

    static sofa::type::Mat<NumberOfNodesInElement, ElementDimension, Real> gradientShapeFunctions(const sofa::type::Vec<ElementDimension, Real>& q)
    {
        return {
            {1. / 4. * (-1 + q[1]), 1. / 4. * (-1 + q[0])},
            {1. / 4. * ( 1 - q[1]), 1. / 4. * (-1 - q[0])},
            {1. / 4. * ( 1 + q[1]), 1. / 4. * ( 1 + q[0])},
            {1. / 4. * (-1 - q[1]), 1. / 4. * ( 1 - q[0])}
        };
    }

    static std::array<QuadraturePointAndWeight, 3> quadraturePoints()
    {
        static sofa::type::Vec<ElementDimension, Real> q0(std::sqrt(2./3.), 0.);
        static sofa::type::Vec<ElementDimension, Real> q1(-1/std::sqrt(6.), -1./std::sqrt(2.));
        static sofa::type::Vec<ElementDimension, Real> q2(-1/std::sqrt(6.), 1./std::sqrt(2.));

        return {
            std::make_pair(q0, 4./3.),
            std::make_pair(q1, 4./3.),
            std::make_pair(q2, 4./3.),
        };
    }
};

}

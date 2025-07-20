#pragma once
#include <Elasticity/Elements.h>

namespace elasticity
{

template<>
inline sofa::type::vector<sofa::topology::Element<sofa::geometry::Edge>> getElementSequence(sofa::core::topology::BaseMeshTopology& topology)
{
    return topology.getEdges();
}

template <class Coord>
struct Volume<sofa::geometry::Edge, Coord>
{
    static SReal compute(const std::array<Coord, sofa::geometry::Edge::NumberOfNodes>& nodesCoordinates)
    {
        return sofa::geometry::Edge::length(nodesCoordinates[0], nodesCoordinates[1]);
    }
};

template<>
constexpr sofa::Size getDimension<sofa::geometry::Edge>()
{
    return 1;
}

template <class DataTypes>
struct ReferenceElement<sofa::geometry::Edge, DataTypes>
{
    using Coord = sofa::Coord_t<DataTypes>;
    using ElementType = sofa::geometry::Edge;
    static constexpr sofa::Size NumberOfNodesInElement = ElementType::NumberOfNodes;

    constexpr static std::array<Coord, NumberOfNodesInElement> nodes = []()
    {
        std::array<Coord, NumberOfNodesInElement> nodes;
        DataTypes::set(nodes[0], -1, 0, 0);
        DataTypes::set(nodes[1], 1, 0, 0);
        return nodes;
    }();
};

}

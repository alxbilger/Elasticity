#pragma once
#include <Elasticity/FiniteElement.h>

namespace elasticity
{

template <class DataTypes>
struct FiniteElement<sofa::geometry::Quad, DataTypes>
{
    using Coord = sofa::Coord_t<DataTypes>;
    using Real = sofa::Real_t<DataTypes>;
    using ElementType = sofa::geometry::Quad;
    using TopologyElement = sofa::topology::Element<ElementType>;
    static constexpr sofa::Size NumberOfNodesInElement = ElementType::NumberOfNodes;
    static constexpr sofa::Size ElementDimension = 2;

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
};

}

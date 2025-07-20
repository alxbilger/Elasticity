#pragma once
#include <Elasticity/Elements.h>

namespace elasticity
{

template<>
inline sofa::type::vector<sofa::topology::Element<sofa::geometry::Quad>> getElementSequence(sofa::core::topology::BaseMeshTopology& topology)
{
    return topology.getQuads();
}

template <class Coord>
struct Volume<sofa::geometry::Quad, Coord>
{
    static SReal compute(const std::array<Coord, sofa::geometry::Quad::NumberOfNodes>& nodesCoordinates)
    {
        return sofa::geometry::Quad::area(nodesCoordinates[0], nodesCoordinates[1], nodesCoordinates[2], nodesCoordinates[3]);
    }
};

template<>
constexpr sofa::Size getDimension<sofa::geometry::Quad>()
{
    return 2;
}

template <class DataTypes>
struct ReferenceElement<sofa::geometry::Quad, DataTypes>
{
    using Coord = sofa::Coord_t<DataTypes>;
    using ElementType = sofa::geometry::Quad;
    static constexpr sofa::Size NumberOfNodesInElement = ElementType::NumberOfNodes;

    constexpr static std::array<Coord, NumberOfNodesInElement> nodes = []()
    {
        std::array<Coord, NumberOfNodesInElement> nodes;
        DataTypes::set(nodes[0], -1, -1, 0);
        DataTypes::set(nodes[1], 1, -1, 0);
        DataTypes::set(nodes[2], 1, 1, 0);
        DataTypes::set(nodes[3], -1, 1, 0);
        return nodes;
    }();
};

}

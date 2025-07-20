#pragma once
#include <Elasticity/Elements.h>

namespace elasticity
{

template<>
inline sofa::type::vector<sofa::topology::Element<sofa::geometry::Hexahedron>> getElementSequence(sofa::core::topology::BaseMeshTopology& topology)
{
    return topology.getHexahedra();
}

template <class Coord>
struct Volume<sofa::geometry::Hexahedron, Coord>
{
    static SReal compute(const std::array<Coord, sofa::geometry::Hexahedron::NumberOfNodes>& nodesCoordinates)
    {
        return sofa::geometry::Hexahedron::volume(nodesCoordinates[0], nodesCoordinates[1], nodesCoordinates[2], nodesCoordinates[3], nodesCoordinates[4], nodesCoordinates[5], nodesCoordinates[6], nodesCoordinates[7]);
    }
};

template<>
constexpr sofa::Size getDimension<sofa::geometry::Hexahedron>()
{
    return 3;
}

template <class DataTypes>
struct ReferenceElement<sofa::geometry::Hexahedron, DataTypes>
{
    using Coord = sofa::Coord_t<DataTypes>;
    using ElementType = sofa::geometry::Hexahedron;
    static constexpr sofa::Size NumberOfNodesInElement = ElementType::NumberOfNodes;

    constexpr static std::array<Coord, NumberOfNodesInElement> nodes = []()
    {
        std::array<Coord, NumberOfNodesInElement> nodes;
        DataTypes::set(nodes[0], -1, -1, -1);
        DataTypes::set(nodes[1], 1, -1, -1);
        DataTypes::set(nodes[2], 1, 1, -1);
        DataTypes::set(nodes[3], -1, 1, -1);
        DataTypes::set(nodes[4], -1, -1, 1);
        DataTypes::set(nodes[5], 1, -1, 1);
        DataTypes::set(nodes[6], 1, 1, 1);
        DataTypes::set(nodes[7], -1, 1, 1);
        return nodes;
    }();
};
}

#pragma once
#include <Elasticity/Elements.h>

namespace elasticity
{

template<>
inline sofa::type::vector<sofa::topology::Element<sofa::geometry::Triangle>> getElementSequence(sofa::core::topology::BaseMeshTopology& topology)
{
    return topology.getTriangles();
}

template <class Coord>
struct Volume<sofa::geometry::Triangle, Coord>
{
    static SReal compute(const std::array<Coord, sofa::geometry::Triangle::NumberOfNodes>& nodesCoordinates)
    {
        return sofa::geometry::Triangle::area(nodesCoordinates[0], nodesCoordinates[1], nodesCoordinates[2]);
    }
};

template<>
constexpr sofa::Size getDimension<sofa::geometry::Triangle>()
{
    return 2;
}

template <class DataTypes>
struct ReferenceElement<sofa::geometry::Triangle, DataTypes>
{
    using Coord = sofa::Coord_t<DataTypes>;
    using ElementType = sofa::geometry::Triangle;
    static constexpr sofa::Size NumberOfNodesInElement = ElementType::NumberOfNodes;

    constexpr static std::array<Coord, NumberOfNodesInElement> nodes = []()
    {
        std::array<Coord, NumberOfNodesInElement> nodes;
        DataTypes::set(nodes[0], 0, 0, 0);
        DataTypes::set(nodes[1], 1, 0, 0);
        DataTypes::set(nodes[2], 0, 1, 0);
        return nodes;
    }();
};

}

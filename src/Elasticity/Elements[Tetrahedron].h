#pragma once
#include <Elasticity/Elements.h>

namespace elasticity
{

template<>
inline sofa::type::vector<sofa::topology::Element<sofa::geometry::Tetrahedron>> getElementSequence(sofa::core::topology::BaseMeshTopology& topology)
{
    return topology.getTetrahedra();
}

template <class Coord>
struct Volume<sofa::geometry::Tetrahedron, Coord>
{
    static SReal compute(const std::array<Coord, sofa::geometry::Tetrahedron::NumberOfNodes>& nodesCoordinates)
    {
        return sofa::geometry::Tetrahedron::volume(nodesCoordinates[0], nodesCoordinates[1], nodesCoordinates[2], nodesCoordinates[3]);
    }
};

template<>
constexpr sofa::Size getDimension<sofa::geometry::Tetrahedron>()
{
    return 3;
}

template <class DataTypes>
struct ReferenceElement<sofa::geometry::Tetrahedron, DataTypes>
{
    using Coord = sofa::Coord_t<DataTypes>;
    using ElementType = sofa::geometry::Tetrahedron;
    static constexpr sofa::Size NumberOfNodesInElement = ElementType::NumberOfNodes;

    constexpr static std::array<Coord, NumberOfNodesInElement> nodes = []()
    {
        std::array<Coord, NumberOfNodesInElement> nodes;
        DataTypes::set(nodes[0], 0, 0, 0);
        DataTypes::set(nodes[1], 1, 0, 0);
        DataTypes::set(nodes[2], 0, 1, 0);
        DataTypes::set(nodes[3], 0, 0, 1);
        return nodes;
    }();
};

}

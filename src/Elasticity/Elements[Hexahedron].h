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
inline sofa::Size getDimension<sofa::geometry::Hexahedron>()
{
    return 3;
}

}

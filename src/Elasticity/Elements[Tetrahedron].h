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
inline sofa::Size getDimension<sofa::geometry::Tetrahedron>()
{
    return 3;
}

}

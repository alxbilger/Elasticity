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
inline sofa::Size getDimension<sofa::geometry::Edge>()
{
    return 1;
}

}

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
inline sofa::Size getDimension<sofa::geometry::Quad>()
{
    return 2;
}

}

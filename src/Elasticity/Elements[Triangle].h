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
inline sofa::Size getDimension<sofa::geometry::Triangle>()
{
    return 2;
}

}

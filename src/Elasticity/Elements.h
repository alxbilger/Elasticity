#pragma once
#include <sofa/core/topology/BaseMeshTopology.h>

namespace elasticity
{

template<class ElementType>
sofa::type::vector<sofa::topology::Element<ElementType>> getElementSequence(sofa::core::topology::BaseMeshTopology& topology);

template <class ElementType, class Coord>
struct Volume;








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

}

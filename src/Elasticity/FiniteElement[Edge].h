#pragma once
#include <Elasticity/FiniteElement.h>

namespace elasticity
{

template <class DataTypes>
struct FiniteElement<sofa::geometry::Edge, DataTypes>
{
    FINITEELEMENT_HEADER(sofa::geometry::Edge, DataTypes, 1);

    constexpr static std::array<Coord, NumberOfNodesInElement> referenceElementNodes = []()
    {
        std::array<Coord, NumberOfNodesInElement> nodes;
        DataTypes::set(nodes[0], -1, 0, 0);
        DataTypes::set(nodes[1], 1, 0, 0);
        return nodes;
    }();

    static sofa::type::vector<TopologyElement> getElementSequence(sofa::core::topology::BaseMeshTopology& topology)
    {
        return topology.getEdges();
    }

    static sofa::type::Mat<NumberOfNodesInElement, ElementDimension, Real> gradientShapeFunctions(const sofa::type::Vec<ElementDimension, Real>& q)
    {
        SOFA_UNUSED(q);
        return {{-1/2}, {1/2}};
    }

    static std::array<QuadraturePointAndWeight, 1> quadraturePoints()
    {
        static sofa::type::Vec<ElementDimension, Real> q0(static_cast<Real>(0));
        return {
            std::make_pair(q0, 2.)
        };
    }

};

}

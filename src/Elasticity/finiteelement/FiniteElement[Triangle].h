#pragma once
#include <Elasticity/finiteelement/FiniteElement.h>

namespace elasticity
{

template <class DataTypes>
struct FiniteElement<sofa::geometry::Triangle, DataTypes>
{
    FINITEELEMENT_HEADER(sofa::geometry::Triangle, DataTypes, 2);
    static_assert(spatial_dimensions > 1, "Triangles cannot be defined in 1D");

    constexpr static std::array<ReferenceCoord, NumberOfNodesInElement> referenceElementNodes {{
        {0, 0},
        {1, 0},
        {0, 1}}};

    static const sofa::type::vector<TopologyElement>& getElementSequence(sofa::core::topology::BaseMeshTopology& topology)
    {
        return topology.getTriangles();
    }

    static constexpr sofa::type::Vec<NumberOfNodesInElement, Real> shapeFunctions(const ReferenceCoord& q)
    {
        // Reference triangle: (x,y) with x>=0, y>=0, x+y<=1
        // Barycentric P1:
        // N0 = 1 - x - y, N1 = x, N2 = y
        return {
            static_cast<Real>(1) - q[0] - q[1],
            q[0],
            q[1]
        };
    }

    static constexpr sofa::type::Mat<NumberOfNodesInElement, TopologicalDimension, Real> gradientShapeFunctions(const sofa::type::Vec<TopologicalDimension, Real>& q)
    {
        SOFA_UNUSED(q);
        return {
            {-1, -1},
            {1, 0},
            {0, 1}
        };
    }

    static constexpr std::array<QuadraturePointAndWeight, 1> quadraturePoints()
    {
        return {
            std::make_pair(sofa::type::Vec<TopologicalDimension, Real>(1./3., 1./3.), 1./2.)
        };
    }
};

}

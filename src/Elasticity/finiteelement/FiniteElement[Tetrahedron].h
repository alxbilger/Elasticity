#pragma once
#include <Elasticity/finiteelement/FiniteElement.h>

namespace elasticity
{

template <class DataTypes>
struct FiniteElement<sofa::geometry::Tetrahedron, DataTypes>
{
    FINITEELEMENT_HEADER(sofa::geometry::Tetrahedron, DataTypes, 3);
    static_assert(spatial_dimensions == 3, "Tetrahedrons are only defined in 3D");

    constexpr static std::array<ReferenceCoord, NumberOfNodesInElement> referenceElementNodes {{
        {0, 0, 0},
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1}
    }};

    static const sofa::type::vector<TopologyElement>& getElementSequence(sofa::core::topology::BaseMeshTopology& topology)
    {
        return topology.getTetrahedra();
    }

    static constexpr sofa::type::Vec<NumberOfNodesInElement, Real> shapeFunctions(const ReferenceCoord& q)
    {
        // Reference tetra: (x,y,z) with x,y,z>=0, x+y+z<=1
        // P1 barycentric:
        // N0 = 1-x-y-z, N1 = x, N2 = y, N3 = z
        return {
            static_cast<Real>(1) - q[0] - q[1] - q[2],
            q[0],
            q[1],
            q[2]
        };
    }

    static constexpr sofa::type::Mat<NumberOfNodesInElement, TopologicalDimension, Real> gradientShapeFunctions(const sofa::type::Vec<TopologicalDimension, Real>& q)
    {
        SOFA_UNUSED(q);
        return {
            {-1, -1, -1},
            {1, 0, 0},
            {0, 1, 0},
            {0, 0, 1}
        };
    }

    static constexpr std::array<QuadraturePointAndWeight, 1> quadraturePoints()
    {
        constexpr sofa::type::Vec<TopologicalDimension, Real> q0(1./4., 1./4., 1./4.);
        constexpr std::array<QuadraturePointAndWeight, 1> q { std::make_pair(q0, 1./6.) };
        return q;
    }
};

}

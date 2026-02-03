#pragma once
#include <Elasticity/finiteelement/FiniteElement.h>

namespace elasticity
{

template <class DataTypes>
struct FiniteElement<sofa::geometry::Hexahedron, DataTypes>
{
    FINITEELEMENT_HEADER(sofa::geometry::Hexahedron, DataTypes, 3);
    static_assert(spatial_dimensions == 3, "Hexahedrons are only defined in 3D");

    constexpr static std::array<ReferenceCoord, NumberOfNodesInElement> referenceElementNodes {{
        {-1, -1, -1},
        {1, -1, -1},
        {1, 1, -1},
        {-1, 1, -1},
        {-1, -1, 1},
        {1, -1, 1},
        {1, 1, 1},
        {-1, 1, 1},
    }};

    static const sofa::type::vector<TopologyElement>& getElementSequence(sofa::core::topology::BaseMeshTopology& topology)
    {
        return topology.getHexahedra();
    }

    static constexpr sofa::type::Vec<NumberOfNodesInElement, Real> shapeFunctions(const ReferenceCoord& q)
    {
        // Reference hexa: (x,y,z) in [-1,1]^3, trilinear (Q1)
        const Real x = q[0];
        const Real y = q[1];
        const Real z = q[2];

        sofa::type::Vec<NumberOfNodesInElement, Real> N(sofa::type::NOINIT);
        for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
        {
            const auto& [xref, yref, zref] = referenceElementNodes[i];
            N[i] = static_cast<Real>(1.0 / 8.0)
                 * (static_cast<Real>(1) + x * xref)
                 * (static_cast<Real>(1) + y * yref)
                 * (static_cast<Real>(1) + z * zref);
        }
        return N;
    }

    static constexpr sofa::type::Mat<NumberOfNodesInElement, TopologicalDimension, Real> gradientShapeFunctions(const sofa::type::Vec<TopologicalDimension, Real>& q)
    {
        const auto [x, y, z] = q;
        sofa::type::Mat<NumberOfNodesInElement, TopologicalDimension, Real> gradient(sofa::type::NOINIT);
        using Line = typename sofa::type::Mat<NumberOfNodesInElement, TopologicalDimension, Real>::Line;

        for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
        {
            const auto& [xref, yref, zref] = referenceElementNodes[i];
            gradient[i] = 1./8. * Line(
                xref * (1 + y * yref) * (1 + z * zref),
                yref * (1 + x * xref) * (1 + z * zref),
                zref * (1 + x * xref) * (1 + y * yref));
        }

        return gradient;
    }

    static constexpr std::array<QuadraturePointAndWeight, 4> quadraturePoints()
    {
        constexpr Real sqrt2_3 = 0.816496580928; //sqrt(2./3.)
        constexpr Real sqrt3 = 1.73205080757; //sqrt(3.)

        constexpr sofa::type::Vec<TopologicalDimension, Real> q0(0., sqrt2_3, -1./sqrt3);
        constexpr sofa::type::Vec<TopologicalDimension, Real> q1(0., -sqrt2_3, -1./sqrt3);
        constexpr sofa::type::Vec<TopologicalDimension, Real> q2(-sqrt2_3, 0., 1./sqrt3);
        constexpr sofa::type::Vec<TopologicalDimension, Real> q3(sqrt2_3, 0., 1./sqrt3);

        constexpr std::array<QuadraturePointAndWeight, 4> q {
            std::make_pair(q0, static_cast<Real>(2)),
            std::make_pair(q1, static_cast<Real>(2)),
            std::make_pair(q2, static_cast<Real>(2)),
            std::make_pair(q3, static_cast<Real>(2))
        };

        return q;
    }
};
}

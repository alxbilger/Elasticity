#pragma once
#include <Elasticity/FiniteElement.h>

namespace elasticity
{

template <class DataTypes>
struct FiniteElement<sofa::geometry::Hexahedron, DataTypes>
{
    FINITEELEMENT_HEADER(sofa::geometry::Hexahedron, DataTypes, 3);
    static_assert(spatial_dimensions == 3, "Hexahedrons are only defined in 3D");

    constexpr static std::array<Coord, NumberOfNodesInElement> referenceElementNodes = []()
    {
        std::array<Coord, NumberOfNodesInElement> nodes;
        DataTypes::set(nodes[0], -1, -1, -1);
        DataTypes::set(nodes[1], 1, -1, -1);
        DataTypes::set(nodes[2], 1, 1, -1);
        DataTypes::set(nodes[3], -1, 1, -1);
        DataTypes::set(nodes[4], -1, -1, 1);
        DataTypes::set(nodes[5], 1, -1, 1);
        DataTypes::set(nodes[6], 1, 1, 1);
        DataTypes::set(nodes[7], -1, 1, 1);
        return nodes;
    }();

    static sofa::type::vector<TopologyElement> getElementSequence(sofa::core::topology::BaseMeshTopology& topology)
    {
        return topology.getHexahedra();
    }

    static sofa::type::Mat<NumberOfNodesInElement, ElementDimension, Real> gradientShapeFunctions(const sofa::type::Vec<ElementDimension, Real>& q)
    {
        const auto [x, y, z] = q;
        sofa::type::Mat<NumberOfNodesInElement, ElementDimension, Real> gradient(sofa::type::NOINIT);
        using Line = typename sofa::type::Mat<NumberOfNodesInElement, ElementDimension, Real>::Line;

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

    static std::array<QuadraturePointAndWeight, 8> quadraturePoints()
    {
        // static sofa::type::Vec<ElementDimension, Real> q0(0., std::sqrt(2./3.), -1./std::sqrt(3.));
        // static sofa::type::Vec<ElementDimension, Real> q1(0., -std::sqrt(2./3.), -1./std::sqrt(3.));
        // static sofa::type::Vec<ElementDimension, Real> q2(-std::sqrt(2./3.), 0., 1./std::sqrt(3.));
        // static sofa::type::Vec<ElementDimension, Real> q3(std::sqrt(2./3.), 0., 1./std::sqrt(3.));
        //
        // static std::array<QuadraturePointAndWeight, 4> q {
        //     std::make_pair(q0, static_cast<Real>(2)),
        //     std::make_pair(q1, static_cast<Real>(2)),
        //     std::make_pair(q2, static_cast<Real>(2)),
        //     std::make_pair(q3, static_cast<Real>(2))
        // };
        //
        // return q;

        static Real k = 1. / std::sqrt(3.);

        static sofa::type::Vec<ElementDimension, Real> q0(k, k, k);
        static sofa::type::Vec<ElementDimension, Real> q1(k, k, -k);
        static sofa::type::Vec<ElementDimension, Real> q2(k, -k, k);
        static sofa::type::Vec<ElementDimension, Real> q3(k, -k, -k);
        static sofa::type::Vec<ElementDimension, Real> q4(-k, k, k);
        static sofa::type::Vec<ElementDimension, Real> q5(-k, k, -k);
        static sofa::type::Vec<ElementDimension, Real> q6(-k, -k, k);
        static sofa::type::Vec<ElementDimension, Real> q7(-k, -k, -k);

        static std::array<QuadraturePointAndWeight, 8> q {
            std::make_pair(q0, static_cast<Real>(1)),
            std::make_pair(q1, static_cast<Real>(1)),
            std::make_pair(q2, static_cast<Real>(1)),
            std::make_pair(q3, static_cast<Real>(1)),
            std::make_pair(q4, static_cast<Real>(1)),
            std::make_pair(q5, static_cast<Real>(1)),
            std::make_pair(q6, static_cast<Real>(1)),
            std::make_pair(q7, static_cast<Real>(1))
        };

        return q;
    }
};
}

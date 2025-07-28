#pragma once
#include <Elasticity/finiteelement/FiniteElement.h>

namespace elasticity
{

template <class DataTypes>
struct FiniteElement<sofa::geometry::Edge, DataTypes>
{
    FINITEELEMENT_HEADER(sofa::geometry::Edge, DataTypes, 1);

    constexpr static std::array<ReferenceCoord, NumberOfNodesInElement> referenceElementNodes {{ReferenceCoord{-1}, ReferenceCoord{1}}};

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

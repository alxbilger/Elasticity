#pragma once
#include <sofa/core/topology/BaseMeshTopology.h>

namespace elasticity
{

template<class ElementType>
sofa::type::vector<sofa::topology::Element<ElementType>> getElementSequence(sofa::core::topology::BaseMeshTopology& topology);

template <class ElementType, class Coord>
struct Volume;

template<class ElementType>
sofa::Size getDimension();

}

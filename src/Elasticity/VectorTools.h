#pragma once

#include <sofa/topology/Element.h>

namespace elasticity
{

template<class ElementType, class VecCoord>
std::array<typename VecCoord::value_type, ElementType::NumberOfNodes>
extractNodesVectorFromGlobalVector(const sofa::topology::Element<ElementType>& element, const VecCoord& vector)
{
    std::array<typename VecCoord::value_type, ElementType::NumberOfNodes> nodes;
    for (sofa::Size i = 0; i < ElementType::NumberOfNodes; ++i)
    {
        nodes[i] = vector[element[i]];
    }
    return nodes;
}

}

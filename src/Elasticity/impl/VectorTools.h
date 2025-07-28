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

template<sofa::Size N, typename ValueType>
ValueType traceFromVoigtTensor(const sofa::type::Vec<N, ValueType>& voigtTensor)
{
    static const sofa::Size TensorSize = (-1 + std::sqrt(1 + 8 * N)) / 2;
    ValueType trace = 0;
    for (sofa::Size i = 0; i < TensorSize; ++i)
    {
        trace += voigtTensor[i];
    }
    return trace;
}

}

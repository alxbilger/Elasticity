#pragma once

#include <sofa/topology/Element.h>

namespace elasticity
{

namespace detail
{

template<class ElementType, class VecCoord, std::size_t... I>
auto extractRefNodesVectorFromGlobalVector(const sofa::topology::Element<ElementType>& element, const VecCoord& vector, std::index_sequence<I...>)
{
    using Coord = typename VecCoord::value_type;
    using CoordRef = std::reference_wrapper<const Coord>;
    return std::array<CoordRef, ElementType::NumberOfNodes>{ std::cref(vector[element[I]])... };
}

template<class ElementType, class VecCoord, std::size_t... I>
auto extractNodesVectorFromGlobalVector(const sofa::topology::Element<ElementType>& element, const VecCoord& vector, std::index_sequence<I...>)
{
    using Coord = typename VecCoord::value_type;
    return std::array<Coord, ElementType::NumberOfNodes>{ vector[element[I]]... };
}

}

template<class ElementType, class VecCoord>
auto extractRefNodesVectorFromGlobalVector(
    const sofa::topology::Element<ElementType>& element, const VecCoord& vector)
{
    return detail::extractRefNodesVectorFromGlobalVector(element, vector, std::make_index_sequence<ElementType::NumberOfNodes>{});
}

template<class ElementType, class VecCoord>
std::array<typename VecCoord::value_type, ElementType::NumberOfNodes>
extractNodesVectorFromGlobalVector(const sofa::topology::Element<ElementType>& element, const VecCoord& vector)
{
    return detail::extractNodesVectorFromGlobalVector(element, vector, std::make_index_sequence<ElementType::NumberOfNodes>{});
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

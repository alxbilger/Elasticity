#pragma once
#include <sofa/config.h>

namespace elasticity::symmetric_tensor
{

/// The number of independent elements in a symmetric 2nd-order tensor of size (N x N)
template<sofa::Size N>
constexpr sofa::Size NumberOfIndependentElements = N * (N + 1) / 2;

}

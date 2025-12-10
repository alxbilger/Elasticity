#pragma once
#include <Elasticity/impl/trait.h>
#include <sofa/helper/OptionsGroup.h>

namespace elasticity
{

constexpr std::string_view parallelComputeStrategy = "parallel";
constexpr std::string_view parallelUnsequencedComputeStrategy = "parallel_unsequenced";
constexpr std::string_view sequencedComputeStrategy = "sequenced";
constexpr std::string_view unsequencedComputeStrategy = "unsequenced";

MAKE_SELECTABLE_ITEMS(ComputeStrategy,
    sofa::helper::Item{parallelComputeStrategy, "The algorithm is executed in parallel"},
    sofa::helper::Item{parallelUnsequencedComputeStrategy, "The algorithm is executed in parallel and may be vectorized"},
    sofa::helper::Item{sequencedComputeStrategy, "The algorithm is executed sequentially"},
    sofa::helper::Item{unsequencedComputeStrategy, "The algorithm may be vectorized"},
);

}

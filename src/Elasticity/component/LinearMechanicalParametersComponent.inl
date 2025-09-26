#pragma once
#include <Elasticity/component/LinearMechanicalParametersComponent.h>

namespace elasticity
{

template <class TReal>
LinearMechanicalParametersComponent<TReal>::LinearMechanicalParametersComponent()
    : d_poissonRatio(initData(&d_poissonRatio, static_cast<TReal>(0.45), "poissonRatio", "Poisson's ratio"))
    , d_youngModulus(initData(&d_youngModulus, static_cast<TReal>(1e6), "youngModulus", "Young's modulus"))
{
}

}

#pragma once

#include <Elasticity/component/material/DirichletMaterial.h>
#include <Elasticity/component/HyperelasticMaterial.inl>

namespace elasticity
{

template <class DataTypes>
auto DirichletMaterial<DataTypes>::firstPiolaKirchhoffStress(const DeformationGradient& F)
-> StressTensor
{
    return static_cast<Real>(2) * F;
}

template <class DataTypes>
auto DirichletMaterial<DataTypes>::materialTangentModulus(const DeformationGradient& F) -> StressJacobian
{
    StressJacobian dPdF;
    return dPdF;
}

}  // namespace elasticity

#pragma once

#include <Elasticity/component/material/DirichletMaterial.h>
#include <Elasticity/impl/ElasticityTensor.h>

#include <Elasticity/component/HyperelasticMaterial.inl>

namespace elasticity
{

template <class DataTypes>
DirichletMaterial<DataTypes>::DirichletMaterial()
{

}

template <class DataTypes>
auto DirichletMaterial<DataTypes>::firstPiolaKirchhoffStress(const DeformationGradient& F)
-> StressTensor
{
    return static_cast<Real>(2) * F;
}

template <class DataTypes>
auto DirichletMaterial<DataTypes>::jacobianFirstPiolaKirchhoffStress() -> StressJacobian
{
    StressJacobian dPdF;
    return dPdF;
}

}  // namespace elasticity

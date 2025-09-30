#pragma once

#include <Elasticity/component/material/NeoHookeanMaterial.h>
#include <Elasticity/impl/MatrixTools.h>

#include <Elasticity/component/HyperelasticMaterial.inl>

namespace elasticity
{

template <class DataTypes>
auto NeoHookeanMaterial<DataTypes>::firstPiolaKirchhoffStress(const DeformationGradient& F)
-> StressTensor
{
    const auto J = elasticity::determinant(F);
    const auto J_inv = 1 / J;

    // derivative of J with respect to F
    const auto dJdF = elasticity::adjugate(F).transposed();

    return m_mu * (F - J_inv *  dJdF) + (m_lambda * log(J) * J_inv) * dJdF;
}

template <class DataTypes>
auto NeoHookeanMaterial<DataTypes>::jacobianFirstPiolaKirchhoffStress(const DeformationGradient& F) -> StressJacobian
{
    // derivative of J with respect to F
    const auto dJdF = elasticity::adjugate(F).transposed();

    const auto dJdF_vec = elasticity::flatten(dJdF);


    StressJacobian dPdF;
    return dPdF;
}

}  // namespace elasticity

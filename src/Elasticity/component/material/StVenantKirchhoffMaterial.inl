#pragma once

#include <Elasticity/component/material/StVenantKirchhoffMaterial.h>
#include <Elasticity/impl/ElasticityTensor.h>

#include <Elasticity/component/HyperelasticMaterial.inl>

namespace elasticity
{

template <class DataTypes>
auto StVenantKirchhoffMaterial<DataTypes>::firstPiolaKirchhoffStress(const DeformationGradient& F)
-> StressTensor
{
    static const auto& I = sofa::type::Mat<spatial_dimensions, spatial_dimensions, Real>::Identity();

    // Right Cauchy-Green deformation tensor
    const auto C = F.transposed() * F;

    // Green-Lagrangian strain tensor
    const auto E = 0.5 * (C - I);

    // Second Piola-Kirchhoff stress tensor
    const auto S = m_lambda * sofa::type::trace(E) * I + 2 * m_mu * E;

    // First Piola-Kirchhoff stress tensor
    const auto P = F * S;

    return P;
}

template <class DataTypes>
auto StVenantKirchhoffMaterial<DataTypes>::jacobianFirstPiolaKirchhoffStress(const DeformationGradient& F) -> StressJacobian
{
    StressJacobian dPdF;
    return dPdF;
}

}  // namespace elasticity

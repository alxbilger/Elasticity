#pragma once

#include <Elasticity/component/material/StVenantKirchhoffMaterial.h>
#include <Elasticity/impl/ElasticityTensor.h>

#include <Elasticity/component/HyperelasticMaterial.inl>

namespace elasticity
{

template <class DataTypes>
auto StVenantKirchhoffMaterial<DataTypes>::secondPiolaKirchhoffStress(const DeformationGradient& C)
-> StressTensor
{
    static const auto& I = sofa::type::Mat<spatial_dimensions, spatial_dimensions, Real>::Identity();

    // Green-Lagrangian strain tensor
    const auto E = 0.5 * (C - I);

    // Second Piola-Kirchhoff stress tensor
    return m_lambda * sofa::type::trace(E) * I + static_cast<Real>(2) * m_mu * E;
}

template <class DataTypes>
auto StVenantKirchhoffMaterial<DataTypes>::elasticityTensor(const DeformationGradient& F) -> StressJacobian
{
    SOFA_UNUSED(F);
    StressJacobian elasticityTensor;

    elasticityTensor.fill(
        [mu = m_mu, lambda = m_lambda](sofa::Index i, sofa::Index j, sofa::Index k, sofa::Index l)
        {
            return mu * (kroneckerDelta(i, k) * kroneckerDelta(j, l) + kroneckerDelta(i, l) * kroneckerDelta(j, k)) +
                        lambda * kroneckerDelta(i, j) * kroneckerDelta(k, l);
        });

    return elasticityTensor;
}

}  // namespace elasticity

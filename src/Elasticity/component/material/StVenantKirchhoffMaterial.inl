#pragma once

#include <Elasticity/component/material/StVenantKirchhoffMaterial.h>
#include <Elasticity/component/PK2HyperelasticMaterial.inl>
#include <sofa/helper/ScopedAdvancedTimer.h>

namespace elasticity
{

template <class DataTypes>
auto StVenantKirchhoffMaterial<DataTypes>::secondPiolaKirchhoffStress(Strain<DataTypes>& strain)
-> StressTensor
{
    static const auto& I = sofa::type::Mat<spatial_dimensions, spatial_dimensions, Real>::Identity();

    // Green-Lagrangian strain tensor
    const auto& E = strain.getGreenLagrangeTensor();

    // Second Piola-Kirchhoff stress tensor
    return m_lambda * sofa::type::trace(E) * I + static_cast<Real>(2) * m_mu * E;
}

template <class DataTypes>
auto StVenantKirchhoffMaterial<DataTypes>::elasticityTensor(Strain<DataTypes>& strain) -> ElasticityTensor
{
    SCOPED_TIMER_TR("elasticityTensor");
    SOFA_UNUSED(strain);

    return makeIsotropicElasticityTensor<DataTypes>(m_mu, m_lambda);
}

}  // namespace elasticity

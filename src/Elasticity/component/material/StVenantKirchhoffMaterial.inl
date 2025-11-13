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

    return ElasticityTensor(
        [mu = m_mu, lambda = m_lambda](sofa::Index i, sofa::Index j, sofa::Index k, sofa::Index l)
        {
            return mu * (kroneckerDelta<Real>(i, k) * kroneckerDelta<Real>(j, l) + kroneckerDelta<Real>(i, l) * kroneckerDelta<Real>(j, k)) +
                        lambda * kroneckerDelta<Real>(i, j) * kroneckerDelta<Real>(k, l);
        });
}

}  // namespace elasticity

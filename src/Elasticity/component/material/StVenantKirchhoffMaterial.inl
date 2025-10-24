#pragma once

#include <Elasticity/component/material/StVenantKirchhoffMaterial.h>
#include <Elasticity/impl/ElasticityTensor.h>

#include <Elasticity/component/HyperelasticMaterial.inl>

namespace elasticity
{

template <class DataTypes>
auto StVenantKirchhoffMaterial<DataTypes>::secondPiolaKirchhoffStress(const DeformationGradient& F)
-> StressTensor
{
    static const auto& I = sofa::type::Mat<spatial_dimensions, spatial_dimensions, Real>::Identity();

    // Right Cauchy-Green deformation tensor
    const auto C = F.transposed() * F;

    // Green-Lagrangian strain tensor
    const auto E = 0.5 * (C - I);

    // Second Piola-Kirchhoff stress tensor
    return m_lambda * sofa::type::trace(E) * I + 2 * m_mu * E;
}

template <class DataTypes>
auto StVenantKirchhoffMaterial<DataTypes>::elasticityTensor(const DeformationGradient& F) -> StressJacobian
{
    StressJacobian C;

    for (std::size_t i = 0; i < spatial_dimensions; ++i)
    {
        for (std::size_t j = 0; j < spatial_dimensions; ++j)
        {
            for (std::size_t k = 0; k < spatial_dimensions; ++k)
            {
                for (std::size_t l = 0; l < spatial_dimensions; ++l)
                {
                    C(i, j, k, l) =
                        m_mu * (kroneckerDelta(i, k) * kroneckerDelta(j, l) + kroneckerDelta(i, l) * kroneckerDelta(j, k))
                        + m_lambda * kroneckerDelta(i, j) * kroneckerDelta(k, l);
                }
            }
        }
    }

    return C;
}

}  // namespace elasticity

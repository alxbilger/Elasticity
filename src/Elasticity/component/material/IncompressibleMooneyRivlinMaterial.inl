#pragma once

#include <Elasticity/component/material/IncompressibleMooneyRivlinMaterial.h>
#include <Elasticity/impl/MatrixTools.h>

#include <Elasticity/component/HyperelasticMaterial.inl>

namespace elasticity
{

template <class DataTypes>
IncompressibleMooneyRivlinMaterial<DataTypes>::IncompressibleMooneyRivlinMaterial()
    : m_mu10(initData(&m_mu10, static_cast<Real>(1e3), "mu10",
                      "Material constant associated to the first invariant"))
    , m_mu01(initData(&m_mu01, static_cast<Real>(1e3), "mu01",
                      "Material constant associated to the second invariant"))
{}

template <class DataTypes>
auto IncompressibleMooneyRivlinMaterial<DataTypes>::secondPiolaKirchhoffStress(const RightCauchyGreenTensor& C) -> StressTensor
{
    static const auto& I = sofa::type::Mat<spatial_dimensions, spatial_dimensions, Real>::Identity();

    const auto trC = sofa::type::trace(C);

    const auto mu10 = m_mu10.getValue();
    const auto mu01 = m_mu01.getValue();

    return 2 * mu10 * I + 2 * mu01 * (trC * I - C);
}

template <class DataTypes>
auto IncompressibleMooneyRivlinMaterial<DataTypes>::elasticityTensor(const RightCauchyGreenTensor& C) -> ElasticityTensor
{
    return ElasticityTensor(
        [mu01 = m_mu01.getValue()](sofa::Index i, sofa::Index j, sofa::Index k, sofa::Index l)
        {
            return 2 * mu01 * (
                2 * kroneckerDelta<Real>(i, j) * kroneckerDelta<Real>(k, l)
                - kroneckerDelta<Real>(i, k) * kroneckerDelta<Real>(j, l)
                - kroneckerDelta<Real>(i, l) * kroneckerDelta<Real>(j, k)
                );
        });
}

}  // namespace elasticity

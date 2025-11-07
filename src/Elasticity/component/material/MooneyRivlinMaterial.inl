#pragma once

#include <Elasticity/component/material/MooneyRivlinMaterial.h>
#include <Elasticity/impl/MatrixTools.h>

#include <Elasticity/component/HyperelasticMaterial.inl>

namespace elasticity
{

template <class DataTypes>
MooneyRivlinMaterial<DataTypes>::MooneyRivlinMaterial()
    : m_mu10(initData(&m_mu10, static_cast<Real>(1e3), "mu10",
                      "Material constant associated to the first invariant"))
    , m_mu01(initData(&m_mu01, static_cast<Real>(1e3), "mu01",
                      "Material constant associated to the second invariant"))
    , m_bulkModulus(initData(&m_bulkModulus, static_cast<Real>(1e6), "bulkModulus", "Bulk modulus"))
{}

template <class DataTypes>
auto MooneyRivlinMaterial<DataTypes>::secondPiolaKirchhoffStress(Strain<DataTypes>& strain) -> StressTensor
{
    static constexpr auto& I = Strain<DataTypes>::identity;
    constexpr auto dim_1 = static_cast<Real>(1) / spatial_dimensions;
    const auto& C = strain.getRightCauchyGreenTensor();
    const auto C_1 = inverse(C);
    const auto J = strain.getDeterminantDeformationGradient();

    const auto invariant1 = strain.getInvariant1();
    const auto invariant2 = strain.getInvariant2();

    const auto mu10 = m_mu10.getValue();
    const auto mu01 = m_mu01.getValue();
    const auto k = m_bulkModulus.getValue();

    return static_cast<Real>(2) * mu10 * pow(J, -static_cast<Real>(2) * dim_1) * (I - dim_1 * invariant1 * C_1)
        + static_cast<Real>(2) * mu01 * pow(J, -static_cast<Real>(4) * dim_1) * (invariant1 * I - C - static_cast<Real>(2) * dim_1 * invariant2 * C_1)
        + (k / static_cast<Real>(2)) * log(J) * C_1
    ;

}

template <class DataTypes>
auto MooneyRivlinMaterial<DataTypes>::elasticityTensor(Strain<DataTypes>& strain) -> ElasticityTensor
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

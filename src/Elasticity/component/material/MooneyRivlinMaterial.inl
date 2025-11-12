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
    , m_bulkModulus(initData(&m_bulkModulus, static_cast<Real>(1e2), "bulkModulus", "Bulk modulus"))
{}

template <class DataTypes>
auto MooneyRivlinMaterial<DataTypes>::secondPiolaKirchhoffStress(Strain<DataTypes>& strain) -> StressTensor
{
    const auto& C = strain.getRightCauchyGreenTensor();
    const auto C_1 = inverse(C);

    const auto J = strain.getDeterminantDeformationGradient();
    assert(J > 0);

    const auto S_isochoric = [this, J, &strain, &C_1, &C]()
    {
        static constexpr auto& I = Strain<DataTypes>::identity;
        static constexpr Real dim_1 = static_cast<Real>(1) / static_cast<Real>(spatial_dimensions);

        const auto mu10 = m_mu10.getValue();
        const auto mu01 = m_mu01.getValue();

        const auto invariant1 = strain.getInvariant1();
        const auto invariant2 = strain.getInvariant2();

        const auto S_mu_10 = pow(J, -static_cast<Real>(2) * dim_1) * (I - dim_1 * invariant1 * C_1);
        const auto S_mu_01 = pow(J, -static_cast<Real>(4) * dim_1) * (invariant1 * I - C - static_cast<Real>(2) * dim_1 * invariant2 * C_1);

        return static_cast<Real>(2) * (mu10 * S_mu_10 + mu01 * S_mu_01);
    }();

    const auto S_volumetric = [this, J, &C_1]()
    {
        const auto bulk = m_bulkModulus.getValue();
        return bulk * log(J) * C_1;
    }();

    return S_isochoric + S_volumetric;
}

template <class DataTypes>
auto MooneyRivlinMaterial<DataTypes>::elasticityTensor(Strain<DataTypes>& strain) -> ElasticityTensor
{
    static constexpr Real dim_1 = static_cast<Real>(1) / static_cast<Real>(spatial_dimensions);
    const auto& C = strain.getRightCauchyGreenTensor();
    const auto J = strain.getDeterminantDeformationGradient();
    const auto logJ = log(J);
    const auto J_2dim = pow(J, -static_cast<Real>(2) * dim_1);
    const auto J_4dim = pow(J, -static_cast<Real>(4) * dim_1);
    const auto C_1 = elasticity::inverse(C);
    const auto I1 = strain.getInvariant1();
    const auto I2 = strain.getInvariant2();

    const auto mu01 = m_mu01.getValue();
    const auto mu10 = m_mu10.getValue();
    const auto bulk = m_bulkModulus.getValue();

    return ElasticityTensor(
        [&](sofa::Index i, sofa::Index j, sofa::Index k, sofa::Index l)
        {
            //derivative of C^{-1} with respect to C
            const auto dC_1dC = -static_cast<Real>(0.5) * (C_1(i, k) * C_1(l, j) + C_1(i, l) * C_1(k, j));

            const Real C_mu_10 =
                -dim_1 * J_2dim * (
                    kroneckerDelta<Real>(i, j) * C_1(k, l) +
                    (kroneckerDelta<Real>(k, l) - dim_1 * C_1(k, l) * I1) * C_1(i, j) +
                    I1 * dC_1dC
                );

            const Real C_mu_01 = J_4dim * (
                -static_cast<Real>(2) * dim_1 * C_1(k, l) * (I1 * kroneckerDelta<Real>(i, j) - C(i, j) - static_cast<Real>(2) * dim_1 * C_1(i, j) * I2)
                + kroneckerDelta<Real>(i, j) * kroneckerDelta<Real>(k, l)
                - kroneckerDelta<Real>(i, k) * kroneckerDelta<Real>(j, l)
                - static_cast<Real>(2) * dim_1 * (dC_1dC * I2 + C_1(i, j) * (I1 * kroneckerDelta<Real>(k, l) - C(k, l)))
            );

            const Real C_isochoric = static_cast<Real>(4) * (mu10 * C_mu_10 + mu01 * C_mu_01);
            const Real C_volumetric = bulk * (C_1(l, k) * C_1(i, j)) + logJ * dC_1dC;

            return C_isochoric + C_volumetric;
        });
}

}  // namespace elasticity

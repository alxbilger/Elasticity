#pragma once

#include <Elasticity/component/material/NeoHookeanMaterial.h>
#include <Elasticity/impl/MatrixTools.h>

#include <Elasticity/component/HyperelasticMaterial.inl>

namespace elasticity
{

template <class DataTypes>
auto NeoHookeanMaterial<DataTypes>::secondPiolaKirchhoffStress(const RightCauchyGreenTensor& C) -> StressTensor
{
    static const auto& I = sofa::type::Mat<spatial_dimensions, spatial_dimensions, Real>::Identity();

    const DeformationGradient C_1 = inverse(C);
    const Real J = determinantSquareMatrix(C);

    return m_mu * (I - C_1) + m_lambda * std::log(J) * C_1;
}

template <class DataTypes>
auto NeoHookeanMaterial<DataTypes>::elasticityTensor(const RightCauchyGreenTensor& C) -> StressJacobian
{
    StressJacobian elasticityTensor;

    const RightCauchyGreenTensor C_1 = elasticity::inverse(C);
    const Real J = elasticity::determinantSquareMatrix(C);
    const Real logJ = std::log(J);

    elasticityTensor.fill(
        [mu = m_mu, lambda = m_lambda, &C_1, logJ](sofa::Index i, sofa::Index j, sofa::Index k, sofa::Index l)
        {
            return (mu - lambda * logJ) * (C_1(i, k) * C_1(l, j) + C_1(i, l) * C_1(k, j))
                + lambda * C_1(l, k) * C_1(i, j);
        });

    return elasticityTensor;
}

}  // namespace elasticity

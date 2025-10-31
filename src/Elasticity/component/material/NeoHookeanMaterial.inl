#pragma once

#include <Elasticity/component/material/NeoHookeanMaterial.h>
#include <Elasticity/impl/MatrixTools.h>

#include <Elasticity/component/HyperelasticMaterial.inl>

namespace elasticity
{

template <class DataTypes>
auto NeoHookeanMaterial<DataTypes>::secondPiolaKirchhoffStress(const DeformationGradient& C) -> StressTensor
{
    static const auto& I = sofa::type::Mat<spatial_dimensions, spatial_dimensions, Real>::Identity();

    const DeformationGradient C_1 = inverse(C);
    const Real J = determinantSquareMatrix(C);

    return m_mu * (I - C_1) + m_lambda * std::log(J) * C_1;
}

template <class DataTypes>
auto NeoHookeanMaterial<DataTypes>::elasticityTensor(const DeformationGradient& C) -> StressJacobian
{
    StressJacobian elasticityTensor;

    const DeformationGradient C_1 = elasticity::inverse(C);
    const Real J = elasticity::determinantSquareMatrix(C);
    const Real logJ = std::log(J);

    for (std::size_t i = 0; i < spatial_dimensions; ++i)
    {
        for (std::size_t j = 0; j < spatial_dimensions; ++j)
        {
            for (std::size_t k = 0; k < spatial_dimensions; ++k)
            {
                for (std::size_t l = 0; l < spatial_dimensions; ++l)
                {
                    elasticityTensor(i, j, k, l) = (m_mu - m_lambda * logJ) * (C_1(i, k) * C_1(l, j) + C_1(i, l) * C_1(k, j))
                        + m_lambda * C_1(l, k) * C_1(i, j);
                }
            }
        }
    }

    return elasticityTensor;
}

}  // namespace elasticity

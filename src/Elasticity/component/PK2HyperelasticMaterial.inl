#pragma once

#include <Elasticity/component/PK2HyperelasticMaterial.h>
#include <Elasticity/component/HyperelasticMaterial.inl>

namespace elasticity
{

template <class TDataTypes>
auto PK2HyperelasticMaterial<TDataTypes>::firstPiolaKirchhoffStress(const DeformationGradient& F) -> StressTensor
{
    const auto C = F.transposed() * F;
    const auto S = secondPiolaKirchhoffStress(C);
    return F * S;
}

template <class TDataTypes>
auto PK2HyperelasticMaterial<TDataTypes>::materialTangentModulus(const DeformationGradient& F) -> StressJacobian
{
    using HyperelasticMaterial<TDataTypes>::spatial_dimensions;
    using HyperelasticMaterial<TDataTypes>::kroneckerDelta;

    StressJacobian A;

    const auto rightCauchyGreenTensor = F.transposed() * F;
    const auto C = elasticityTensor(rightCauchyGreenTensor);
    const auto S = secondPiolaKirchhoffStress(rightCauchyGreenTensor);

    A.fill([&F, &C, &S](sofa::Index i, sofa::Index j, sofa::Index k, sofa::Index l)
    {
        auto A_ijkl = kroneckerDelta(i,k) * S(l, j);
        for (std::size_t q = 0; q < spatial_dimensions; ++q)
        {
            for (std::size_t r = 0; r < spatial_dimensions; ++r)
            {
                A_ijkl += F(i, q) * C(q, j, l, r) * F(k, r);
            }
        }
        return A_ijkl;
    });
    return A;
}

}  // namespace elasticity

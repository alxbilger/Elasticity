#pragma once
#include <Elasticity/component/HyperelasticMaterial.h>
#include <Elasticity/impl/MatrixTools.h>

namespace elasticity
{

template <class DataTypes>
void HyperelasticMaterial<DataTypes>::init()
{
    BaseObject::init();

    if (!this->isComponentStateInvalid())
    {
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Valid);
    }
}

template <class DataTypes>
auto HyperelasticMaterial<DataTypes>::firstPiolaKirchhoffStress(const DeformationGradient& F) -> StressTensor
{
    const auto S = secondPiolaKirchhoffStress(F);
    return F * S;
}

template <class DataTypes>
auto HyperelasticMaterial<DataTypes>::materialTangentModulus(const DeformationGradient& F) -> StressJacobian
{
    StressJacobian A;
    const auto C = elasticityTensor(F);
    const auto S = secondPiolaKirchhoffStress(F);

    for (std::size_t i = 0; i < spatial_dimensions; ++i)
    {
        for (std::size_t j = 0; j < spatial_dimensions; ++j)
        {
            for (std::size_t k = 0; k < spatial_dimensions; ++k)
            {
                for (std::size_t l = 0; l < spatial_dimensions; ++l)
                {
                    auto& A_ijkl = A(i, j, k, l);
                    A_ijkl = kroneckerDelta(i,k) * S(l, j);
                    for (std::size_t q = 0; q < spatial_dimensions; ++q)
                    {
                        for (std::size_t r = 0; r < spatial_dimensions; ++r)
                        {
                            A_ijkl += 2 * F(i, q) * C(q, j, l, r) * F(k, r);
                        }
                    }
                }
            }
        }
    }
    return A;
}

template <class DataTypes>
auto HyperelasticMaterial<DataTypes>::invariant1(const DeformationGradient& F) -> Real
{
    return elasticity::squaredFrobeniusNorm(F);
}

template <class DataTypes>
auto HyperelasticMaterial<DataTypes>::invariant2(const DeformationGradient& F)  -> Real
{
    const auto C = F.transposed() * F;
    const auto trC2 = elasticity::squaredFrobeniusNorm(C);
    return 0.5 * (invariant1(F) - trC2);
}

template <class DataTypes>
auto HyperelasticMaterial<DataTypes>::invariant3(const DeformationGradient& F)  -> Real
{
    return std::pow(elasticity::determinantSquareMatrix(F), static_cast<Real>(2));
}

}  // namespace elasticity

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

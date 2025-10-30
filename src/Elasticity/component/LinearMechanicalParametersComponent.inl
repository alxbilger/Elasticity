#pragma once
#include <Elasticity/component/LinearMechanicalParametersComponent.h>
#include <Elasticity/impl/ElasticityTensor.h>

namespace elasticity
{

template <class DataTypes>
LinearMechanicalParametersComponent<DataTypes>::LinearMechanicalParametersComponent()
: d_poissonRatio(initData(&d_poissonRatio, static_cast<Real>(0.45), "poissonRatio",
    "Poisson's ratio: represents the material's ability to undergo deformation in directions orthogonal to the applied stress"))
, d_youngModulus(initData(&d_youngModulus, static_cast<Real>(1e6), "youngModulus",
    "Young's modulus: describes the material's stiffness"))
{
    this->addUpdateCallback("toLameCoefficients", {&this->d_youngModulus, &this->d_poissonRatio},
    [this](const sofa::core::DataTracker& )
    {
        setLameCoefficients();
        return this->getComponentState();
    }, {});
    setLameCoefficients();
}

template <class DataTypes>
void LinearMechanicalParametersComponent<DataTypes>::setLameCoefficients()
{
    std::tie(m_lambda, m_mu) = elasticity::toLameParameters<DataTypes>(
        this->d_youngModulus.getValue(), this->d_poissonRatio.getValue());
}

}  // namespace elasticity

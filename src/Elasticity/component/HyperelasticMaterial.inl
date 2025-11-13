#pragma once
#include <Elasticity/component/HyperelasticMaterial.h>

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

}  // namespace elasticity

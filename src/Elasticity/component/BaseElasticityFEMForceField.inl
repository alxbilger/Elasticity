#pragma once

#include <elasticity/component/BaseElasticityFEMForceField.h>
#include <sofa/core/ObjectFactory.h>

namespace elasticity
{

template <class DataTypes>
BaseElasticityFEMForceField<DataTypes>::BaseElasticityFEMForceField()
    : l_topology(initLink("topology", "Link to a topology containing elements"))
{
}

template <class DataTypes>
void BaseElasticityFEMForceField<DataTypes>::init()
{
    sofa::core::behavior::ForceField<DataTypes>::init();

    if (!this->isComponentStateInvalid())
    {
        this->validateTopology();
    }

    if (!this->isComponentStateInvalid())
    {
        selectFEMTypes();
    }
}

template <class DataTypes>
void BaseElasticityFEMForceField<DataTypes>::validateTopology()
{
    if (l_topology.empty())
    {
        msg_info() << "Link to Topology container should be set to ensure right behavior. First "
                      "Topology found in current context will be used.";
        l_topology.set(this->getContext()->getMeshTopologyLink());
    }

    if (l_topology == nullptr)
    {
        msg_error() << "No topology component found at path: " << this->l_topology.getLinkedPath()
                    << ", nor in current context: " << this->getContext()->name
                    << ". Object must have a BaseMeshTopology. "
                    << "The list of available BaseMeshTopology components is: "
                    << sofa::core::ObjectFactory::getInstance()
                           ->listClassesDerivedFrom<sofa::core::topology::BaseMeshTopology>();
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
    }
}

}  // namespace elasticity

#pragma once
#include <Elasticity/BaseLinearSmallStrainFEMForceField.h>
#include <sofa/core/ObjectFactory.h>

namespace elasticity
{

template <class DataTypes>
BaseLinearSmallStrainFEMForceField<DataTypes>::BaseLinearSmallStrainFEMForceField()
    : l_topology(initLink("topology", "Link to a topology containing elements")),
      d_poissonRatio(
          initData(&d_poissonRatio, static_cast<Real>(0.45), "poissonRatio", "Poisson's ratio")),
      d_youngModulus(
          initData(&d_youngModulus, static_cast<Real>(1e6), "youngModulus", "Young's modulus"))
{
    static std::string groupName = "Mechanical parameter";
    d_poissonRatio.setGroup(groupName);
    d_youngModulus.setGroup(groupName);
}

template <class DataTypes>
void BaseLinearSmallStrainFEMForceField<DataTypes>::init()
{
    Inherit1::init();

    if (!this->isComponentStateInvalid())
    {
        validateTopology();
    }

    if (!this->isComponentStateInvalid())
    {
        precomputeElementStiffness();
    }

    if (!this->isComponentStateInvalid())
    {
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Valid);
    }
}

template <class DataTypes>
void BaseLinearSmallStrainFEMForceField<DataTypes>::validateTopology()
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

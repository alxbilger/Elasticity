#include <Elasticity/component/TopologyAccessor.h>
#include <sofa/core/ObjectFactory.h>

namespace elasticity
{


TopologyAccessor::TopologyAccessor()
    : l_topology(initLink("topology", "Link to a topology containing elements"))
{
}

void TopologyAccessor::init()
{
    sofa::core::objectmodel::BaseObject::init();

    if (!this->isComponentStateInvalid())
    {
        this->validateTopology();
    }
}

void TopologyAccessor::validateTopology()
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

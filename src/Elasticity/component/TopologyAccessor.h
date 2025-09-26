#pragma once

#include <Elasticity/config.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/topology/BaseMeshTopology.h>

namespace elasticity
{

class ELASTICITY_API TopologyAccessor : public virtual sofa::core::objectmodel::BaseObject
{
public:
    SOFA_CLASS(TopologyAccessor, sofa::core::objectmodel::BaseObject);

    void init() override;

    /// The topology will give access to the elements
    sofa::SingleLink<TopologyAccessor, sofa::core::topology::BaseMeshTopology,
        sofa::BaseLink::FLAG_STOREPATH | sofa::BaseLink::FLAG_STRONGLINK> l_topology;

protected:

    TopologyAccessor();

    /**
     * Ensure a link to a valid topology. Without a topology, this component cannot have access
     * to the list of elements.
     */
    void validateTopology();
};

}

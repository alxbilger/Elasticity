#pragma once

#include <Elasticity/config.h>
#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/topology/BaseMeshTopology.h>

namespace elasticity
{

template <class DataTypes>
class BaseElasticityFEMForceField : public sofa::core::behavior::ForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(BaseElasticityFEMForceField, DataTypes),
               SOFA_TEMPLATE(sofa::core::behavior::ForceField, DataTypes));

    /// The topology will give access to the elements
    sofa::SingleLink<BaseElasticityFEMForceField, sofa::core::topology::BaseMeshTopology,
        sofa::BaseLink::FLAG_STOREPATH | sofa::BaseLink::FLAG_STRONGLINK> l_topology;

    void init() override;

protected:
    BaseElasticityFEMForceField();

    /**
     * Ensure a link to a valid topology. Without a topology, this force field cannot have access
     * to the list of elements.
     */
    void validateTopology();

    virtual void selectFEMTypes() = 0;
};

#if !defined(ELASTICITY_COMPONENT_BASE_ELASTICITY_FEM_FORCE_FIELD_CPP)
extern template class ELASTICITY_API BaseElasticityFEMForceField<sofa::defaulttype::Vec1Types>;
extern template class ELASTICITY_API BaseElasticityFEMForceField<sofa::defaulttype::Vec2Types>;
extern template class ELASTICITY_API BaseElasticityFEMForceField<sofa::defaulttype::Vec3Types>;
#endif

}  // namespace elasticity

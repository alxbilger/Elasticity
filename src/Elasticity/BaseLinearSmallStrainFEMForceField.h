#pragma once

#include <Elasticity/config.h>
#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/topology/BaseMeshTopology.h>

namespace elasticity
{

template <class DataTypes>
class BaseLinearSmallStrainFEMForceField : public sofa::core::behavior::ForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(BaseLinearSmallStrainFEMForceField, DataTypes),
               SOFA_TEMPLATE(sofa::core::behavior::ForceField, DataTypes));

private:
    using DataVecCoord = sofa::DataVecDeriv_t<DataTypes>;
    using DataVecDeriv = sofa::DataVecDeriv_t<DataTypes>;
    using Real = sofa::Real_t<DataTypes>;

public:
    void init() override;

    /// The topology will give access to the elements
    sofa::SingleLink<BaseLinearSmallStrainFEMForceField, sofa::core::topology::BaseMeshTopology,
                     sofa::BaseLink::FLAG_STOREPATH | sofa::BaseLink::FLAG_STRONGLINK> l_topology;

    sofa::Data<Real> d_poissonRatio;
    sofa::Data<Real> d_youngModulus;

private:
    /**
     * Ensure a link to a valid topology. Without a topology, this force field cannot have access
     * to the list of elements.
     */
    void validateTopology();

protected:

    BaseLinearSmallStrainFEMForceField();

    /**
     * With linear small strain, the element stiffness matrix is constant, so it can be precomputed.
     */
    virtual void precomputeElementStiffness() = 0;
};

#if !defined(ELASTICITY_BASELINEARSMALLSTRAINFEMFORCEFIELD_CPP)
extern template class ELASTICITY_API BaseLinearSmallStrainFEMForceField<sofa::defaulttype::Vec1Types>;
extern template class ELASTICITY_API BaseLinearSmallStrainFEMForceField<sofa::defaulttype::Vec2Types>;
extern template class ELASTICITY_API BaseLinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types>;
#endif
}  // namespace elasticity

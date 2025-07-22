#pragma once
#include <Elasticity/LinearSmallStrainFEMForceField.h>
#include <Elasticity/BaseLinearSmallStrainFEMForceField.inl>

namespace elasticity
{

template <class DataTypes>
void LinearSmallStrainFEMForceField<DataTypes>::addForce(
    const sofa::core::MechanicalParams* mparams, DataVecDeriv& f, const DataVecCoord& x,
    const DataVecDeriv& v)
{
}

template <class DataTypes>
void LinearSmallStrainFEMForceField<DataTypes>::addDForce(
    const sofa::core::MechanicalParams* mparams, DataVecDeriv& df, const DataVecDeriv& dx)
{
}

template <class DataTypes>
void LinearSmallStrainFEMForceField<DataTypes>::buildStiffnessMatrix(
    sofa::core::behavior::StiffnessMatrix* matrix)
{
}

template <class DataTypes>
SReal LinearSmallStrainFEMForceField<DataTypes>::getPotentialEnergy(
    const sofa::core::MechanicalParams*, const DataVecCoord& x) const
{
    return 0;
}

template <class DataTypes>
void LinearSmallStrainFEMForceField<DataTypes>::precomputeElementStiffness()
{
}

}  // namespace elasticity

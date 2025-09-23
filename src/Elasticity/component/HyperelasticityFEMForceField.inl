#pragma once

#include <Elasticity/component/HyperelasticityFEMForceField.h>
#include <Elasticity/component/BaseElasticityFEMForceField.inl>

namespace elasticity
{

template <class DataTypes>
void HyperelasticityFEMForceField<DataTypes>::addForce(const sofa::core::MechanicalParams* mparams,
                                                       DataVecDeriv& f, const DataVecCoord& x,
                                                       const DataVecDeriv& v)
{
}

template <class DataTypes>
void HyperelasticityFEMForceField<DataTypes>::addDForce(const sofa::core::MechanicalParams* mparams,
                                                        DataVecDeriv& df, const DataVecDeriv& dx)
{
}

template <class DataTypes>
void HyperelasticityFEMForceField<DataTypes>::buildStiffnessMatrix(
    sofa::core::behavior::StiffnessMatrix* matrix)
{
}

template <class DataTypes>
SReal HyperelasticityFEMForceField<DataTypes>::getPotentialEnergy(
    const sofa::core::MechanicalParams*, const DataVecCoord& x) const
{
    return 0;
}

template <class DataTypes>
void HyperelasticityFEMForceField<DataTypes>::selectFEMTypes()
{
}

}

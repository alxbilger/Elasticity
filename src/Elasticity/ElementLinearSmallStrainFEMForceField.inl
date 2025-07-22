#pragma once
#include <Elasticity/ElementLinearSmallStrainFEMForceField.h>
#include <Elasticity/BaseLinearSmallStrainFEMForceField.inl>

#include <Elasticity/LinearFEM.inl>

namespace elasticity
{

template <class DataTypes, class ElementType>
void ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::precomputeElementStiffness()
{
    m_finiteElement.setTopology(l_topology.get());

    auto restPositionAccessor = this->mstate->readRestPositions();
    m_finiteElement.precomputeElementStiffness(restPositionAccessor.ref(), d_youngModulus.getValue(), d_poissonRatio.getValue());
}

template <class DataTypes, class ElementType>
void ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::addForce(
    const sofa::core::MechanicalParams* mparams,
    DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& v)
{
    SOFA_UNUSED(mparams);
    SOFA_UNUSED(v);

    auto forceAccessor = sofa::helper::getWriteOnlyAccessor(f);
    auto positionAccessor = sofa::helper::getReadAccessor(x);
    auto restPositionAccessor = this->mstate->readRestPositions();

    m_finiteElement.addForce(forceAccessor.wref(), positionAccessor.ref(), restPositionAccessor.ref());
}

template <class DataTypes, class ElementType>
void ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::addDForce(
    const sofa::core::MechanicalParams* mparams, DataVecDeriv& df, const DataVecDeriv& dx)
{
    auto dfAccessor = sofa::helper::getWriteAccessor(df);
    auto dxAccessor = sofa::helper::getReadAccessor(dx);
    dfAccessor.resize(dxAccessor.size());

    const Real kFactor = (Real)sofa::core::mechanicalparams::kFactorIncludingRayleighDamping(
        mparams, this->rayleighStiffness.getValue());

    m_finiteElement.addDForce(dfAccessor.wref(), dxAccessor.ref(), kFactor);
}

template <class DataTypes, class ElementType>
void ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::buildStiffnessMatrix(
    sofa::core::behavior::StiffnessMatrix* matrix)
{
    auto dfdx = matrix->getForceDerivativeIn(this->mstate)
                   .withRespectToPositionsIn(this->mstate);

    m_finiteElement.buildStiffnessMatrix(dfdx);
}

template <class DataTypes, class ElementType>
SReal ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::getPotentialEnergy(
    const sofa::core::MechanicalParams*, const DataVecCoord& x) const
{
    return 0;
}

}  // namespace elasticity

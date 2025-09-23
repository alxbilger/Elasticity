#pragma once
#include <Elasticity/impl/FEM/NonLinearFEM.h>

namespace elasticity
{

template <class DataTypes, class ElementType>
NonLinearFEM<DataTypes, ElementType>::NonLinearFEM(sofa::core::topology::BaseMeshTopology* topology)
    : m_topology(topology)
{
}

template <class DataTypes, class ElementType>
void NonLinearFEM<DataTypes, ElementType>::addForce(VecDeriv& force, const VecCoord& position,
                                                    const VecCoord& restPosition)
{
}

template <class DataTypes, class ElementType>
void NonLinearFEM<DataTypes, ElementType>::addDForce(VecDeriv& df, const VecDeriv& dx, Real kFactor) const
{}

template <class DataTypes, class ElementType>
void NonLinearFEM<DataTypes, ElementType>::buildStiffnessMatrix(
    sofa::core::behavior::StiffnessMatrix::Derivative& dfdx) const
{
}

}  // namespace elasticity

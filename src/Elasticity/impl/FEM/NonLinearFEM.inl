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
    if (m_topology == nullptr) return;

    const auto& elements = FiniteElement::getElementSequence(*m_topology);

    for (sofa::Size i = 0; i < elements.size(); ++i)
    {
    }
}

template <class DataTypes, class ElementType>
void NonLinearFEM<DataTypes, ElementType>::addDForce(VecDeriv& df, const VecDeriv& dx, Real kFactor) const
{
    if (m_topology == nullptr) return;

    const auto& elements = FiniteElement::getElementSequence(*m_topology);

    for (sofa::Size i = 0; i < elements.size(); ++i)
    {
    }
}

template <class DataTypes, class ElementType>
void NonLinearFEM<DataTypes, ElementType>::buildStiffnessMatrix(
    sofa::core::behavior::StiffnessMatrix::Derivative& dfdx) const
{
    if (m_topology == nullptr) return;

    const auto& elements = FiniteElement::getElementSequence(*m_topology);

    for (sofa::Size i = 0; i < elements.size(); ++i)
    {
    }
}

}  // namespace elasticity

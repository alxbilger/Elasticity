#pragma once
#include <Elasticity/component/ElementLinearSmallStrainFEMForceField.h>
#include <Elasticity/component/LinearSmallStrainFEMForceField.inl>

#include <Elasticity/impl/FEM/LinearFEM.inl>

namespace elasticity
{

template <class DataTypes, class ElementType>
void ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::selectFEMTypes()
{
    this->template addLinearFEMType<ElementType>();
}

}  // namespace elasticity

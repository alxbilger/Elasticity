#pragma once
#include <Elasticity/ElementLinearSmallStrainFEMForceField.h>
#include <Elasticity/LinearSmallStrainFEMForceField.inl>

#include <Elasticity/LinearFEM.inl>

namespace elasticity
{

template <class DataTypes, class ElementType>
void ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::selectFEMTypes()
{
    this->addLinearFEMType<ElementType>();
}

}  // namespace elasticity

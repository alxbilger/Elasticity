#pragma once
#include <Elasticity/component/prefab/LinearSmallStrainFEMForceField.h>
#include <Elasticity/component/ElementLinearSmallStrainFEMForceField.inl>
#include <Elasticity/component/ElementPrefab.inl>

namespace elasticity
{

template <class DataTypes>
LinearSmallStrainFEMForceField<DataTypes>::LinearSmallStrainFEMForceField()
: d_computeForceStrategy(initData(&d_computeForceStrategy, "computeForceStrategy", std::string("The compute strategy used to compute the forces.\n" + ComputeStrategy::dataDescription()).c_str()))
, d_computeForceDerivStrategy(initData(&d_computeForceDerivStrategy, "computeForceDerivStrategy", std::string("The compute strategy used to compute the forces derivatives.\n" + ComputeStrategy::dataDescription()).c_str()))
, d_elementSpace(initData(&d_elementSpace, static_cast<sofa::Real_t<DataTypes>>(0.125), "elementSpace", "When rendering, the space between elements"))
{
    d_elementSpace.setGroup("Visualization");
}

}  // namespace elasticity

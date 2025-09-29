#pragma once

#include <Elasticity/component/ElementPrefab.h>
#include <Elasticity/component/LinearMechanicalParametersComponent.h>
#include <Elasticity/component/ElementHyperelasticityFEMForceField.h>

#if !defined(ELASTICITY_COMPONENT_PREFAB_HYPERELASTICITY_FEM_FORCEFIELD_CPP)
#include <sofa/defaulttype/VecTypes.h>
#endif

namespace elasticity
{

template <class DataTypes>
using HyperelasticityPrefabParent =
    ElementPrefab<ElementPrefabTrait<ElementHyperelasticityFEMForceField, DataTypes>>;

/**
 * An intermediate class that instantiates other element-specific hyperelasticity components based on
 * the types of elements found in the linked topology.
 *
 * The instantiated components are of type ElementHyperelasticityFEMForceField. Since this class
 * and ElementHyperelasticityFEMForceField both share the same data, they will be linked.
 */
template <class DataTypes>
class HyperelasticityFEMForceField :
    public HyperelasticityPrefabParent<DataTypes>
{
public:
    SOFA_CLASS(HyperelasticityFEMForceField<DataTypes>, HyperelasticityPrefabParent<DataTypes>);
};

#if !defined(ELASTICITY_COMPONENT_PREFAB_HYPERELASTICITY_FEM_FORCEFIELD_CPP)
extern template class ELASTICITY_API HyperelasticityFEMForceField<sofa::defaulttype::Vec1Types>;
extern template class ELASTICITY_API HyperelasticityFEMForceField<sofa::defaulttype::Vec2Types>;
extern template class ELASTICITY_API HyperelasticityFEMForceField<sofa::defaulttype::Vec3Types>;
#endif

}  // namespace elasticity

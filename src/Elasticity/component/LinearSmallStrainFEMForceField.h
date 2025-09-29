#pragma once

#include <Elasticity/component/ElementPrefab.h>
#include <Elasticity/component/LinearMechanicalParametersComponent.h>
#include <Elasticity/component/ElementLinearSmallStrainFEMForceField.h>

#if !defined(ELASTICITY_COMPONENT_LINEAR_SMALL_STRAIN_FEM_FORCEFIELD_CPP)
#include <sofa/defaulttype/VecTypes.h>
#endif

namespace elasticity
{

template <class DataTypes>
using LinearPrefabParent =
    ElementPrefab<ElementPrefabTrait<ElementLinearSmallStrainFEMForceField, DataTypes>>;

/**
 * An intermediate class that instantiates other element-specific hyperelasticity components based on
 * the types of elements found in the linked topology.
 *
 * The instantiated components are of type ElementLinearSmallStrainFEMForceField. Since this class
 * and ElementLinearSmallStrainFEMForceField both share the same data, they will be linked.
 */
template <class DataTypes>
class LinearSmallStrainFEMForceField :
    public LinearPrefabParent<DataTypes>,
    public LinearMechanicalParametersComponent<sofa::Real_t<DataTypes>>
{
public:
    SOFA_CLASS2(
        LinearSmallStrainFEMForceField<DataTypes>,
            LinearPrefabParent<DataTypes>,
            LinearMechanicalParametersComponent<sofa::Real_t<DataTypes>>);
};

#if !defined(ELASTICITY_COMPONENT_LINEAR_SMALL_STRAIN_FEM_FORCEFIELD_CPP)
extern template class ELASTICITY_API LinearSmallStrainFEMForceField<sofa::defaulttype::Vec1Types>;
extern template class ELASTICITY_API LinearSmallStrainFEMForceField<sofa::defaulttype::Vec2Types>;
extern template class ELASTICITY_API LinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types>;
#endif

}  // namespace elasticity

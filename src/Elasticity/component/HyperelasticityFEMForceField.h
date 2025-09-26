#pragma once

#include <Elasticity/component/ElementPrefab.h>
#include <Elasticity/component/LinearMechanicalParametersComponent.h>
#include <Elasticity/component/ElementHyperelasticityFEMForceField.h>

#if !defined(ELASTICITY_COMPONENT_LINEAR_SMALL_STRAIN_FEM_FORCEFIELD_CPP)
#include <sofa/defaulttype/VecTypes.h>
#endif

namespace elasticity
{

template <class DataTypes>
class HyperelasticityFEMForceField :
    public ElementPrefab<ElementPrefabTrait<ElementHyperelasticityFEMForceField, DataTypes>>,
    public LinearMechanicalParametersComponent<sofa::Real_t<DataTypes>>
{
public:
    SOFA_CLASS2(
        HyperelasticityFEMForceField<DataTypes>,
            ElementPrefab<SOFA_TEMPLATE2(ElementPrefabTrait, ElementHyperelasticityFEMForceField, DataTypes)>,
            LinearMechanicalParametersComponent<sofa::Real_t<DataTypes>>);
};

#if !defined(ELASTICITY_COMPONENT_HYPERELASTICITY_FEM_FORCEFIELD_CPP)
extern template class ELASTICITY_API HyperelasticityFEMForceField<sofa::defaulttype::Vec1Types>;
extern template class ELASTICITY_API HyperelasticityFEMForceField<sofa::defaulttype::Vec2Types>;
extern template class ELASTICITY_API HyperelasticityFEMForceField<sofa::defaulttype::Vec3Types>;
#endif

}  // namespace elasticity

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
class LinearSmallStrainFEMForceField :
    public ElementPrefab<ElementPrefabTrait<ElementLinearSmallStrainFEMForceField, DataTypes>>,
    public LinearMechanicalParametersComponent<sofa::Real_t<DataTypes>>
{
public:
    SOFA_CLASS2(
        LinearSmallStrainFEMForceField<DataTypes>,
            ElementPrefab<SOFA_TEMPLATE2(ElementPrefabTrait, ElementLinearSmallStrainFEMForceField, DataTypes)>,
            LinearMechanicalParametersComponent<sofa::Real_t<DataTypes>>);
};

#if !defined(ELASTICITY_COMPONENT_LINEAR_SMALL_STRAIN_FEM_FORCEFIELD_CPP)
extern template class ELASTICITY_API LinearSmallStrainFEMForceField<sofa::defaulttype::Vec1Types>;
extern template class ELASTICITY_API LinearSmallStrainFEMForceField<sofa::defaulttype::Vec2Types>;
extern template class ELASTICITY_API LinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types>;
#endif

}  // namespace elasticity

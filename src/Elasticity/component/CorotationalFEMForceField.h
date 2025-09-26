#pragma once

#include <Elasticity/component/ElementPrefab.h>
#include <Elasticity/component/LinearMechanicalParametersComponent.h>
#include <Elasticity/component/ElementCorotationalFEMForceField.h>

#if !defined(ELASTICITY_COMPONENT_COROTATIONAL_FEM_FORCEFIELD_CPP)
#include <sofa/defaulttype/VecTypes.h>
#endif

namespace elasticity
{

template <class DataTypes>
class CorotationalFEMForceField :
    public ElementPrefab<ElementPrefabTrait<ElementCorotationalFEMForceField, DataTypes>>,
    public LinearMechanicalParametersComponent<sofa::Real_t<DataTypes>>
{
public:
    SOFA_CLASS2(
        CorotationalFEMForceField<DataTypes>,
            ElementPrefab<SOFA_TEMPLATE2(ElementPrefabTrait, ElementCorotationalFEMForceField, DataTypes)>,
            LinearMechanicalParametersComponent<sofa::Real_t<DataTypes>>);
};

#if !defined(ELASTICITY_COMPONENT_COROTATIONAL_FEM_FORCEFIELD_CPP)
// extern template class ELASTICITY_API CorotationalFEMForceField<sofa::defaulttype::Vec1Types>;
extern template class ELASTICITY_API CorotationalFEMForceField<sofa::defaulttype::Vec2Types>;
extern template class ELASTICITY_API CorotationalFEMForceField<sofa::defaulttype::Vec3Types>;
#endif

}  // namespace elasticity

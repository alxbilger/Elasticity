#pragma once

#include <Elasticity/component/ElementPrefab.h>
#include <Elasticity/component/LinearMechanicalParametersComponent.h>
#include <Elasticity/component/ElementCorotationalFEMForceField.h>

#if !defined(ELASTICITY_COMPONENT_PREFAB_COROTATIONAL_FEM_FORCEFIELD_CPP)
#include <sofa/defaulttype/VecTypes.h>
#endif

namespace elasticity
{

template <class DataTypes>
using CorotationalPrefabParent =
    ElementPrefab<ElementPrefabTrait<ElementCorotationalFEMForceField, DataTypes>>;

/**
 * An intermediate class that instantiates other element-specific corotational components based on
 * the types of elements found in the linked topology.
 *
 * The instantiated components are of type ElementCorotationalFEMForceField. Since this class
 * and ElementCorotationalFEMForceField both share the same data, they will be linked.
 */
template <class DataTypes>
class CorotationalFEMForceField :
    public CorotationalPrefabParent<DataTypes>,
    public LinearMechanicalParametersComponent<DataTypes>
{
public:
    SOFA_CLASS2(
        CorotationalFEMForceField<DataTypes>,
            CorotationalPrefabParent<DataTypes>,
            LinearMechanicalParametersComponent<DataTypes>);
};

#if !defined(ELASTICITY_COMPONENT_PREFAB_COROTATIONAL_FEM_FORCEFIELD_CPP)
// extern template class ELASTICITY_API CorotationalFEMForceField<sofa::defaulttype::Vec1Types>;
extern template class ELASTICITY_API CorotationalFEMForceField<sofa::defaulttype::Vec2Types>;
extern template class ELASTICITY_API CorotationalFEMForceField<sofa::defaulttype::Vec3Types>;
#endif

}  // namespace elasticity

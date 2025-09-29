#define ELASTICITY_COMPONENT_PREFAB_HYPERELASTICITY_FEM_FORCEFIELD_CPP

#include <Elasticity/component/prefab/HyperelasticityFEMForceField.inl>
#include <sofa/core/ObjectFactory.h>

namespace elasticity
{

void registerHyperelasticityFEMForceField(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData("Hyperelasticity on linear elements")
        .add< HyperelasticityFEMForceField<sofa::defaulttype::Vec1Types> >()
        .add< HyperelasticityFEMForceField<sofa::defaulttype::Vec2Types> >()
        .add< HyperelasticityFEMForceField<sofa::defaulttype::Vec3Types> >(true));
}

template class ELASTICITY_API HyperelasticityFEMForceField<sofa::defaulttype::Vec1Types>;
template class ELASTICITY_API HyperelasticityFEMForceField<sofa::defaulttype::Vec2Types>;
template class ELASTICITY_API HyperelasticityFEMForceField<sofa::defaulttype::Vec3Types>;

}

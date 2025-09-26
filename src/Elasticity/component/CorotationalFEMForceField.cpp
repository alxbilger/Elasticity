#define ELASTICITY_COMPONENT_COROTATIONAL_FEM_FORCEFIELD_CPP

#include <Elasticity/component/CorotationalFEMForceField.inl>
#include <sofa/core/ObjectFactory.h>

namespace elasticity
{

void registerCorotationalFEMForceField(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData("Hooke's law on linear elements using corotational approach")
        // .add< CorotationalFEMForceField<sofa::defaulttype::Vec1Types> >()
        .add< CorotationalFEMForceField<sofa::defaulttype::Vec2Types> >()
        .add< CorotationalFEMForceField<sofa::defaulttype::Vec3Types> >(true));
}

// template class ELASTICITY_API CorotationalFEMForceField<sofa::defaulttype::Vec1Types>;
template class ELASTICITY_API CorotationalFEMForceField<sofa::defaulttype::Vec2Types>;
template class ELASTICITY_API CorotationalFEMForceField<sofa::defaulttype::Vec3Types>;

}

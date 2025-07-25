#include <Elasticity/CorotationalFEMForceField.inl>

namespace elasticity
{

void registerCorotationalFEMForceField(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData("Hooke's law on linear elements using the corotational approach")
        .add< CorotationalFEMForceField<sofa::defaulttype::Vec3Types> >(true)
        .add< CorotationalFEMForceField<sofa::defaulttype::Vec2Types> >());
}

template class ELASTICITY_API CorotationalFEMForceField<sofa::defaulttype::Vec3Types>;
template class ELASTICITY_API CorotationalFEMForceField<sofa::defaulttype::Vec2Types>;

}

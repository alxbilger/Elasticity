#include <Elasticity/component/HyperelasticityFEMForceField.inl>

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

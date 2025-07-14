#include <Elasticity/TET4LinearSmallStrainFEMForceField.inl>

namespace elasticity
{

void registerTET4LinearSmallStrainFEMForceField(sofa::core::ObjectFactory* factory)
{
  factory->registerObjects(sofa::core::ObjectRegistrationData("Hooke's law on linear tetrahedra assuming small strain")
      .add< TET4LinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types> >(true));
}

template class ELASTICITY_API TET4LinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types>;

}

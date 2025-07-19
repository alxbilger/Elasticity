#include <Elasticity/LinearSmallStrainFEMForceField.inl>

namespace elasticity
{

void registerTET4LinearSmallStrainFEMForceField(sofa::core::ObjectFactory* factory)
{
  factory->registerObjects(sofa::core::ObjectRegistrationData("Hooke's law on linear tetrahedra assuming small strain")
      .add< LinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Tetrahedron> >(true));
}

template class ELASTICITY_API LinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Tetrahedron>;

}

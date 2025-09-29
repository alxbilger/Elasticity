#define ELASTICITY_COMPONENT_PREFAB_LINEAR_SMALL_STRAIN_FEM_FORCEFIELD_CPP

#include <Elasticity/component/prefab/LinearSmallStrainFEMForceField.inl>
#include <sofa/core/ObjectFactory.h>

namespace elasticity
{

void registerLinearSmallStrainFEMForceField(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData("Hooke's law on linear elements assuming small strain")
        .add< LinearSmallStrainFEMForceField<sofa::defaulttype::Vec1Types> >()
        .add< LinearSmallStrainFEMForceField<sofa::defaulttype::Vec2Types> >()
        .add< LinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types> >(true));
}

template class ELASTICITY_API LinearSmallStrainFEMForceField<sofa::defaulttype::Vec1Types>;
template class ELASTICITY_API LinearSmallStrainFEMForceField<sofa::defaulttype::Vec2Types>;
template class ELASTICITY_API LinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types>;

}

#define ELASTICITY_COMPONENT_MATERIAL_INCOMPRESSIBLEMOONEYRIVLIN_CPP

#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/VecTypes.h>

#include <Elasticity/component/material/IncompressibleMooneyRivlinMaterial.inl>

namespace elasticity
{

void registerIncompressibleMooneyRivlinMaterial(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData("Mooney-Rivlin material")
        .add< IncompressibleMooneyRivlinMaterial<sofa::defaulttype::Vec1Types> >()
        .add< IncompressibleMooneyRivlinMaterial<sofa::defaulttype::Vec2Types> >()
        .add< IncompressibleMooneyRivlinMaterial<sofa::defaulttype::Vec3Types> >(true));
}

template class ELASTICITY_API IncompressibleMooneyRivlinMaterial<sofa::defaulttype::Vec1Types>;
template class ELASTICITY_API IncompressibleMooneyRivlinMaterial<sofa::defaulttype::Vec2Types>;
template class ELASTICITY_API IncompressibleMooneyRivlinMaterial<sofa::defaulttype::Vec3Types>;

}

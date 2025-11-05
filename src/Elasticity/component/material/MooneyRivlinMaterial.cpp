#define ELASTICITY_COMPONENT_MATERIAL_MOONEYRIVLIN_CPP

#include <Elasticity/component/material/MooneyRivlinMaterial.inl>

#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/VecTypes.h>

namespace elasticity
{

void registerMooneyRivlinMaterial(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData("Mooney-Rivlin material")
        .add< MooneyRivlinMaterial<sofa::defaulttype::Vec1Types> >()
        .add< MooneyRivlinMaterial<sofa::defaulttype::Vec2Types> >()
        .add< MooneyRivlinMaterial<sofa::defaulttype::Vec3Types> >(true));
}

template class ELASTICITY_API MooneyRivlinMaterial<sofa::defaulttype::Vec1Types>;
template class ELASTICITY_API MooneyRivlinMaterial<sofa::defaulttype::Vec2Types>;
template class ELASTICITY_API MooneyRivlinMaterial<sofa::defaulttype::Vec3Types>;

}

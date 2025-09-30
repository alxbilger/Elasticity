#define ELASTICITY_COMPONENT_MATERIAL_NEOHOOKEANMATERIAL_CPP

#include <Elasticity/component/material/NeoHookeanMaterial.inl>

#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/VecTypes.h>

namespace elasticity
{

void registerNeoHookeanMaterial(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData("Neo-Hookean material")
        .add< NeoHookeanMaterial<sofa::defaulttype::Vec1Types> >()
        .add< NeoHookeanMaterial<sofa::defaulttype::Vec2Types> >()
        .add< NeoHookeanMaterial<sofa::defaulttype::Vec3Types> >(true));
}

template class ELASTICITY_API NeoHookeanMaterial<sofa::defaulttype::Vec1Types>;
template class ELASTICITY_API NeoHookeanMaterial<sofa::defaulttype::Vec2Types>;
template class ELASTICITY_API NeoHookeanMaterial<sofa::defaulttype::Vec3Types>;

}

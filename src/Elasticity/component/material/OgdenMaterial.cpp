#define ELASTICITY_COMPONENT_MATERIAL_OGDEN_CPP

#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/VecTypes.h>

#include <Elasticity/component/material/OgdenMaterial.inl>

namespace elasticity
{

void registerOgdenMaterial(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData("Ogden material")
        .add< OgdenMaterial<sofa::defaulttype::Vec1Types> >()
        .add< OgdenMaterial<sofa::defaulttype::Vec2Types> >()
        .add< OgdenMaterial<sofa::defaulttype::Vec3Types> >(true));
}

template class ELASTICITY_API OgdenMaterial<sofa::defaulttype::Vec1Types>;
template class ELASTICITY_API OgdenMaterial<sofa::defaulttype::Vec2Types>;
template class ELASTICITY_API OgdenMaterial<sofa::defaulttype::Vec3Types>;

}

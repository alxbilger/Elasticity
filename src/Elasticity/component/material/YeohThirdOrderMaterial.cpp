#define ELASTICITY_COMPONENT_MATERIAL_YEOHTHIRDORDERMATERIAL_CPP

#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/VecTypes.h>

#include <Elasticity/component/material/YeohThirdOrderMaterial.inl>

namespace elasticity
{

void registerYeohThirdOrderMaterial(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData("3rd-order Yeoh material")
        .add< YeohThirdOrderMaterial<sofa::defaulttype::Vec1Types> >()
        .add< YeohThirdOrderMaterial<sofa::defaulttype::Vec2Types> >()
        .add< YeohThirdOrderMaterial<sofa::defaulttype::Vec3Types> >(true));
}

template class ELASTICITY_API YeohThirdOrderMaterial<sofa::defaulttype::Vec1Types>;
template class ELASTICITY_API YeohThirdOrderMaterial<sofa::defaulttype::Vec2Types>;
template class ELASTICITY_API YeohThirdOrderMaterial<sofa::defaulttype::Vec3Types>;

}

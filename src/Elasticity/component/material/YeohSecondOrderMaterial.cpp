#define ELASTICITY_COMPONENT_MATERIAL_YEOHSECONDORDERMATERIAL_CPP

#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/VecTypes.h>

#include <Elasticity/component/material/YeohSecondOrderMaterial.inl>

namespace elasticity
{

void registerYeohSecondOrderMaterial(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData("2nd-order Yeoh material")
        .add< YeohSecondOrderMaterial<sofa::defaulttype::Vec1Types> >()
        .add< YeohSecondOrderMaterial<sofa::defaulttype::Vec2Types> >()
        .add< YeohSecondOrderMaterial<sofa::defaulttype::Vec3Types> >(true));
}

template class ELASTICITY_API YeohSecondOrderMaterial<sofa::defaulttype::Vec1Types>;
template class ELASTICITY_API YeohSecondOrderMaterial<sofa::defaulttype::Vec2Types>;
template class ELASTICITY_API YeohSecondOrderMaterial<sofa::defaulttype::Vec3Types>;

}

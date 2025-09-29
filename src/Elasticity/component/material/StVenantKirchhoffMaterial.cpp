#define ELASTICITY_COMPONENT_MATERIAL_ST_VENANT_KIRCHHOFF_MATERIAL_CPP

#include <Elasticity/component/material/StVenantKirchhoffMaterial.inl>

#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/VecTypes.h>

namespace elasticity
{

void registerStVenantKirchhoffMaterial(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData("Saint Venant-Kirchhoff material")
        .add< StVenantKirchhoffMaterial<sofa::defaulttype::Vec1Types> >()
        .add< StVenantKirchhoffMaterial<sofa::defaulttype::Vec2Types> >()
        .add< StVenantKirchhoffMaterial<sofa::defaulttype::Vec3Types> >(true));
}

template class ELASTICITY_API StVenantKirchhoffMaterial<sofa::defaulttype::Vec1Types>;
template class ELASTICITY_API StVenantKirchhoffMaterial<sofa::defaulttype::Vec2Types>;
template class ELASTICITY_API StVenantKirchhoffMaterial<sofa::defaulttype::Vec3Types>;

}

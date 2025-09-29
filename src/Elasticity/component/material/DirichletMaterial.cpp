#define ELASTICITY_COMPONENT_MATERIAL_DIRICHLETMATERIAL_CPP

#include <Elasticity/component/material/DirichletMaterial.inl>

#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/VecTypes.h>

namespace elasticity
{

void registerDirichletMaterial(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData("Dirichlet material")
        .add< DirichletMaterial<sofa::defaulttype::Vec1Types> >()
        .add< DirichletMaterial<sofa::defaulttype::Vec2Types> >()
        .add< DirichletMaterial<sofa::defaulttype::Vec3Types> >(true));
}

template class ELASTICITY_API DirichletMaterial<sofa::defaulttype::Vec1Types>;
template class ELASTICITY_API DirichletMaterial<sofa::defaulttype::Vec2Types>;
template class ELASTICITY_API DirichletMaterial<sofa::defaulttype::Vec3Types>;

}

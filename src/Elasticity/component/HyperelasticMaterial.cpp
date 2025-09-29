#define ELASTICITY_COMPONENT_HYPERELASTIC_MATERIAL_CPP

#include <Elasticity/component/HyperelasticMaterial.inl>
#include <sofa/defaulttype/VecTypes.h>

namespace elasticity
{

template class ELASTICITY_API HyperelasticMaterial<sofa::defaulttype::Vec1Types>;
template class ELASTICITY_API HyperelasticMaterial<sofa::defaulttype::Vec2Types>;
template class ELASTICITY_API HyperelasticMaterial<sofa::defaulttype::Vec3Types>;

}

#define ELASTICITY_COMPONENT_PK2HYPERELASTIC_MATERIAL_CPP
#include <Elasticity/component/PK2HyperelasticMaterial.inl>

namespace elasticity
{

template class ELASTICITY_API PK2HyperelasticMaterial<sofa::defaulttype::Vec1Types>;
template class ELASTICITY_API PK2HyperelasticMaterial<sofa::defaulttype::Vec2Types>;
template class ELASTICITY_API PK2HyperelasticMaterial<sofa::defaulttype::Vec3Types>;

}

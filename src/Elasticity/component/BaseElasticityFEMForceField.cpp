#define ELASTICITY_COMPONENT_BASE_ELASTICITY_FEM_FORCE_FIELD_CPP
#include <Elasticity/component/BaseElasticityFEMForceField.inl>

namespace elasticity
{

template class ELASTICITY_API BaseElasticityFEMForceField<sofa::defaulttype::Vec1Types>;
template class ELASTICITY_API BaseElasticityFEMForceField<sofa::defaulttype::Vec2Types>;
template class ELASTICITY_API BaseElasticityFEMForceField<sofa::defaulttype::Vec3Types>;

}

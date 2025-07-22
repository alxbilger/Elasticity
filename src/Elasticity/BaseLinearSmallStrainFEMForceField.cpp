#define ELASTICITY_BASELINEARSMALLSTRAINFEMFORCEFIELD_CPP
#include <Elasticity/BaseLinearSmallStrainFEMForceField.inl>

namespace elasticity
{

template class ELASTICITY_API BaseLinearSmallStrainFEMForceField<sofa::defaulttype::Vec1Types>;
template class ELASTICITY_API BaseLinearSmallStrainFEMForceField<sofa::defaulttype::Vec2Types>;
template class ELASTICITY_API BaseLinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types>;

}

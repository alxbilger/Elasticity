#define ELASTICITY_COMPONENT_LINEAR_MECHANICAL_PARAMETERS_CPP
#include <Elasticity/component/LinearMechanicalParametersComponent.inl>

namespace elasticity
{

template class ELASTICITY_API LinearMechanicalParametersComponent<sofa::defaulttype::Vec1Types>;
template class ELASTICITY_API LinearMechanicalParametersComponent<sofa::defaulttype::Vec2Types>;
template class ELASTICITY_API LinearMechanicalParametersComponent<sofa::defaulttype::Vec3Types>;

}

#pragma once

#include <Elasticity/config.h>
#include <sofa/core/objectmodel/BaseObject.h>

#if !defined(ELASTICITY_COMPONENT_LINEAR_MECHANICAL_PARAMETERS_CPP)
#include <sofa/defaulttype/VecTypes.h>
#endif

namespace elasticity
{

template<class DataTypes>
class LinearMechanicalParametersComponent : public virtual sofa::core::objectmodel::BaseObject
{
    using Real = sofa::Real_t<DataTypes>;

public:
    SOFA_CLASS(LinearMechanicalParametersComponent<DataTypes>, sofa::core::objectmodel::BaseObject);

    sofa::Data<Real> d_poissonRatio;
    sofa::Data<Real> d_youngModulus;

protected:
    LinearMechanicalParametersComponent();

    // Lamé's coefficients
    Real m_lambda, m_mu;
};

#if !defined(ELASTICITY_COMPONENT_LINEAR_MECHANICAL_PARAMETERS_CPP)
extern template class ELASTICITY_API LinearMechanicalParametersComponent<sofa::defaulttype::Vec1Types>;
extern template class ELASTICITY_API LinearMechanicalParametersComponent<sofa::defaulttype::Vec2Types>;
extern template class ELASTICITY_API LinearMechanicalParametersComponent<sofa::defaulttype::Vec3Types>;
#endif
}  // namespace elasticity

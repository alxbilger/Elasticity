#pragma once

#include <Elasticity/config.h>
#include <sofa/core/objectmodel/BaseObject.h>

namespace elasticity
{

template<class TReal>
class LinearMechanicalParametersComponent : public virtual sofa::core::objectmodel::BaseObject
{
public:
    SOFA_CLASS(LinearMechanicalParametersComponent<TReal>, sofa::core::objectmodel::BaseObject);

    sofa::Data<TReal> d_poissonRatio;
    sofa::Data<TReal> d_youngModulus;

protected:
    LinearMechanicalParametersComponent();
};

#if !defined(ELASTICITY_COMPONENT_LINEAR_MECHANICAL_PARAMETERS_CPP)
extern template class ELASTICITY_API LinearMechanicalParametersComponent<SReal>;
#endif
}  // namespace elasticity

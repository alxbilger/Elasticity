#pragma once

#include <Elasticity/config.h>
#include <Elasticity/component/HyperelasticMaterial.h>
#include <Elasticity/component/LinearMechanicalParametersComponent.h>

#if !defined(ELASTICITY_COMPONENT_MATERIAL_NEOHOOKEANMATERIAL_CPP)
#include <sofa/defaulttype/VecTypes.h>
#endif

namespace elasticity
{

/**
 * Formulation from Bonet, J. and R. D. Wood (2008). Nonlinear continuum mechanics for finite
 * element analysis. Cambridge university press.
 */
template <class DataTypes>
class NeoHookeanMaterial :
    public HyperelasticMaterial<DataTypes>,
    public LinearMechanicalParametersComponent<DataTypes>
{
public:
    SOFA_CLASS2(NeoHookeanMaterial, HyperelasticMaterial<DataTypes>, LinearMechanicalParametersComponent<DataTypes>);

private:
    using Real = sofa::Real_t<DataTypes>;

    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;

    using DeformationGradient = HyperelasticMaterial<DataTypes>::DeformationGradient;
    using StressTensor = HyperelasticMaterial<DataTypes>::StressTensor;
    using StressJacobian = HyperelasticMaterial<DataTypes>::StressJacobian;

    using LinearMechanicalParametersComponent<DataTypes>::m_lambda;
    using LinearMechanicalParametersComponent<DataTypes>::m_mu;

public:
    StressTensor firstPiolaKirchhoffStress(const DeformationGradient& F) override;

    StressJacobian materialTangentModulus(const DeformationGradient& F) override;
};

#if !defined(ELASTICITY_COMPONENT_MATERIAL_NEOHOOKEANMATERIAL_CPP)
extern template class ELASTICITY_API NeoHookeanMaterial<sofa::defaulttype::Vec1Types>;
extern template class ELASTICITY_API NeoHookeanMaterial<sofa::defaulttype::Vec2Types>;
extern template class ELASTICITY_API NeoHookeanMaterial<sofa::defaulttype::Vec3Types>;
#endif
}

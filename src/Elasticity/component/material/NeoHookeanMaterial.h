#pragma once

#include <Elasticity/config.h>
#include <Elasticity/component/PK2HyperelasticMaterial.h>
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
    public PK2HyperelasticMaterial<DataTypes>,
    public LinearMechanicalParametersComponent<DataTypes>
{
public:
    SOFA_CLASS2(NeoHookeanMaterial, PK2HyperelasticMaterial<DataTypes>, LinearMechanicalParametersComponent<DataTypes>);

private:
    using Real = sofa::Real_t<DataTypes>;

    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;

    using DeformationGradient = PK2HyperelasticMaterial<DataTypes>::DeformationGradient;
    using RightCauchyGreenTensor = PK2HyperelasticMaterial<DataTypes>::RightCauchyGreenTensor;
    using StressTensor = PK2HyperelasticMaterial<DataTypes>::StressTensor;
    using ElasticityTensor = PK2HyperelasticMaterial<DataTypes>::ElasticityTensor;

    using LinearMechanicalParametersComponent<DataTypes>::m_lambda;
    using LinearMechanicalParametersComponent<DataTypes>::m_mu;

public:
    StressTensor secondPiolaKirchhoffStress(const RightCauchyGreenTensor& C) override;

    ElasticityTensor elasticityTensor(const RightCauchyGreenTensor& C) override;
};


#if !defined(ELASTICITY_COMPONENT_MATERIAL_NEOHOOKEANMATERIAL_CPP)
extern template class ELASTICITY_API NeoHookeanMaterial<sofa::defaulttype::Vec1Types>;
extern template class ELASTICITY_API NeoHookeanMaterial<sofa::defaulttype::Vec2Types>;
extern template class ELASTICITY_API NeoHookeanMaterial<sofa::defaulttype::Vec3Types>;
#endif
}

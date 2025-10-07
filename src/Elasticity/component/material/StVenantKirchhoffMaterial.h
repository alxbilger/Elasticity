#pragma once

#include <Elasticity/config.h>
#include <Elasticity/component/HyperelasticMaterial.h>
#include <Elasticity/component/LinearMechanicalParametersComponent.h>

#if !defined(ELASTICITY_COMPONENT_MATERIAL_STVENANTKIRCHHOFFMATERIAL_CPP)
#include <sofa/defaulttype/VecTypes.h>
#endif

namespace elasticity
{

template<class DataTypes>
class StVenantKirchhoffMaterial :
    public HyperelasticMaterial<DataTypes>,
    public LinearMechanicalParametersComponent<DataTypes>
{
public:
    SOFA_CLASS2(StVenantKirchhoffMaterial, HyperelasticMaterial<DataTypes>, LinearMechanicalParametersComponent<DataTypes>);

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

    StressJacobian jacobianFirstPiolaKirchhoffStress(const DeformationGradient& F) override;
};

#if !defined(ELASTICITY_COMPONENT_MATERIAL_STVENANTKIRCHHOFFMATERIAL_CPP)
extern template class ELASTICITY_API StVenantKirchhoffMaterial<sofa::defaulttype::Vec1Types>;
extern template class ELASTICITY_API StVenantKirchhoffMaterial<sofa::defaulttype::Vec2Types>;
extern template class ELASTICITY_API StVenantKirchhoffMaterial<sofa::defaulttype::Vec3Types>;
#endif
}

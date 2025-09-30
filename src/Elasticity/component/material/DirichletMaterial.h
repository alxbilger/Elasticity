#pragma once

#include <Elasticity/config.h>
#include <Elasticity/component/HyperelasticMaterial.h>
#include <Elasticity/component/LinearMechanicalParametersComponent.h>

#if !defined(ELASTICITY_COMPONENT_MATERIAL_DIRICHLETMATERIAL_CPP)
#include <sofa/defaulttype/VecTypes.h>
#endif

namespace elasticity
{

template<class DataTypes>
class DirichletMaterial :
    public HyperelasticMaterial<DataTypes>
{
public:
    SOFA_CLASS2(DirichletMaterial, HyperelasticMaterial<DataTypes>, LinearMechanicalParametersComponent<DataTypes>);

private:
    using Real = sofa::Real_t<DataTypes>;

    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;

    using HyperelasticMaterial<DataTypes>::DeformationGradient;
    using HyperelasticMaterial<DataTypes>::StressTensor;
    using HyperelasticMaterial<DataTypes>::StressJacobian;

public:
    StressTensor firstPiolaKirchhoffStress(const DeformationGradient& F) override;

    StressJacobian jacobianFirstPiolaKirchhoffStress(const DeformationGradient& F) override;
};

#if !defined(ELASTICITY_COMPONENT_MATERIAL_DIRICHLETMATERIAL_CPP)
template class ELASTICITY_API DirichletMaterial<sofa::defaulttype::Vec1Types>;
template class ELASTICITY_API DirichletMaterial<sofa::defaulttype::Vec2Types>;
template class ELASTICITY_API DirichletMaterial<sofa::defaulttype::Vec3Types>;
#endif
}

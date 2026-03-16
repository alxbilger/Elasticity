#pragma once

#include <Elasticity/config.h>
#include <Elasticity/component/PK2HyperelasticMaterial.h>

#if !defined(ELASTICITY_COMPONENT_MATERIAL_YeohSecondOrderMaterial_CPP)
#include <sofa/defaulttype/VecTypes.h>
#endif

namespace elasticity
{

template <class DataTypes>
class YeohSecondOrderMaterial :
    public PK2HyperelasticMaterial<DataTypes>
{
public:
    SOFA_CLASS(YeohSecondOrderMaterial, PK2HyperelasticMaterial<DataTypes>);

private:
    using Real = sofa::Real_t<DataTypes>;

    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;

    using DeformationGradient = PK2HyperelasticMaterial<DataTypes>::DeformationGradient;
    using RightCauchyGreenTensor = PK2HyperelasticMaterial<DataTypes>::RightCauchyGreenTensor;
    using StressTensor = PK2HyperelasticMaterial<DataTypes>::StressTensor;
    using ElasticityTensor = PK2HyperelasticMaterial<DataTypes>::ElasticityTensor;

public:
    sofa::Data<Real> m_c1;
    sofa::Data<Real> m_c2;
    sofa::Data<Real> m_bulkModulus1;
    sofa::Data<Real> m_bulkModulus2;

    StressTensor secondPiolaKirchhoffStress(Strain<DataTypes>& strain) override;

    ElasticityTensor elasticityTensor(Strain<DataTypes>& strain) override;

protected:
    YeohSecondOrderMaterial();
};


#if !defined(ELASTICITY_COMPONENT_MATERIAL_YEOHSECONDORDERMATERIAL_CPP)
extern template class ELASTICITY_API YeohSecondOrderMaterial<sofa::defaulttype::Vec1Types>;
extern template class ELASTICITY_API YeohSecondOrderMaterial<sofa::defaulttype::Vec2Types>;
extern template class ELASTICITY_API YeohSecondOrderMaterial<sofa::defaulttype::Vec3Types>;
#endif
}

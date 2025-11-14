#pragma once

#include <Elasticity/config.h>
#include <Elasticity/component/PK2HyperelasticMaterial.h>

#if !defined(ELASTICITY_COMPONENT_MATERIAL_OGDEN_CPP)
#include <sofa/defaulttype/VecTypes.h>
#endif

namespace elasticity
{

template <class DataTypes>
class OgdenMaterial :
    public PK2HyperelasticMaterial<DataTypes>
{
public:
    SOFA_CLASS(OgdenMaterial, PK2HyperelasticMaterial<DataTypes>);

private:
    using Real = sofa::Real_t<DataTypes>;

    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;

    using DeformationGradient = PK2HyperelasticMaterial<DataTypes>::DeformationGradient;
    using RightCauchyGreenTensor = PK2HyperelasticMaterial<DataTypes>::RightCauchyGreenTensor;
    using StressTensor = PK2HyperelasticMaterial<DataTypes>::StressTensor;
    using ElasticityTensor = PK2HyperelasticMaterial<DataTypes>::ElasticityTensor;

public:
    // Ogden material constants
    sofa::Data<Real> m_mu;
    sofa::Data<Real> m_alpha;
    sofa::Data<Real> m_kappa;

    StressTensor secondPiolaKirchhoffStress(Strain<DataTypes>& strain) override;

    ElasticityTensor elasticityTensor(Strain<DataTypes>& strain) override;

protected:
    OgdenMaterial();
};


#if !defined(ELASTICITY_COMPONENT_MATERIAL_OGDEN_CPP)
extern template class ELASTICITY_API OgdenMaterial<sofa::defaulttype::Vec1Types>;
extern template class ELASTICITY_API OgdenMaterial<sofa::defaulttype::Vec2Types>;
extern template class ELASTICITY_API OgdenMaterial<sofa::defaulttype::Vec3Types>;
#endif
}

#pragma once
#include <Elasticity/component/HyperelasticMaterial.h>

namespace elasticity
{

/**
 * A hyperelastic material defined by its second Piola-Kirchhoff stress tensosr and its Lagrangian
 * elasticity tensor.
 */
template <class TDataTypes>
class PK2HyperelasticMaterial : public HyperelasticMaterial<TDataTypes>
{
public:
    SOFA_CLASS(PK2HyperelasticMaterial<TDataTypes>, HyperelasticMaterial<TDataTypes>);
    using DataTypes = TDataTypes;

protected:
    using HyperelasticMaterial<TDataTypes>::DeformationGradient;
    using HyperelasticMaterial<TDataTypes>::RightCauchyGreenTensor;
    using HyperelasticMaterial<TDataTypes>::StressTensor;
    using HyperelasticMaterial<TDataTypes>::ElasticityTensor;
    using HyperelasticMaterial<TDataTypes>::TangentModulus;
    using HyperelasticMaterial<TDataTypes>::spatial_dimensions;
    using HyperelasticMaterial<TDataTypes>::kroneckerDelta;

public:
    StressTensor firstPiolaKirchhoffStress(Strain<DataTypes>& strain) final;
    TangentModulus materialTangentModulus(Strain<DataTypes>& strain) final;

protected:
    virtual StressTensor secondPiolaKirchhoffStress(Strain<DataTypes>& strain) = 0;
    virtual ElasticityTensor elasticityTensor(Strain<DataTypes>& strain) = 0;
};

#if !defined(ELASTICITY_COMPONENT_HYPERELASTIC_MATERIAL_CPP)
extern template class ELASTICITY_API PK2HyperelasticMaterial<sofa::defaulttype::Vec1Types>;
extern template class ELASTICITY_API PK2HyperelasticMaterial<sofa::defaulttype::Vec2Types>;
extern template class ELASTICITY_API PK2HyperelasticMaterial<sofa::defaulttype::Vec3Types>;
#endif

}  // namespace elasticity

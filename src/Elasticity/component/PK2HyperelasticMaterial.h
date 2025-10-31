#pragma once
#include <Elasticity/component/HyperelasticMaterial.h>

namespace elasticity
{

template<class TDataTypes>
class PK2HyperelasticMaterial : public HyperelasticMaterial<TDataTypes>
{
public:
    SOFA_CLASS(PK2HyperelasticMaterial<TDataTypes>, HyperelasticMaterial<TDataTypes>);
    using DataTypes = TDataTypes;

protected:
    using HyperelasticMaterial<TDataTypes>::DeformationGradient;
    using HyperelasticMaterial<TDataTypes>::StressTensor;
    using HyperelasticMaterial<TDataTypes>::StressJacobian;
    using HyperelasticMaterial<TDataTypes>::spatial_dimensions;
    using HyperelasticMaterial<TDataTypes>::kroneckerDelta;

public:
    StressTensor firstPiolaKirchhoffStress(const DeformationGradient& F) final;
    StressJacobian materialTangentModulus(const DeformationGradient& F) final;

protected:
    virtual StressTensor secondPiolaKirchhoffStress(const DeformationGradient& C) = 0;
    virtual StressJacobian elasticityTensor(const DeformationGradient& C) = 0;
};

#if !defined(ELASTICITY_COMPONENT_HYPERELASTIC_MATERIAL_CPP)
extern template class ELASTICITY_API PK2HyperelasticMaterial<sofa::defaulttype::Vec1Types>;
extern template class ELASTICITY_API PK2HyperelasticMaterial<sofa::defaulttype::Vec2Types>;
extern template class ELASTICITY_API PK2HyperelasticMaterial<sofa::defaulttype::Vec3Types>;
#endif

}  // namespace elasticity

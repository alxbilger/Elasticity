#pragma once

#include <Elasticity/config.h>
#include <Elasticity/impl/FullySymmetric4Tensor.h>
#include <Elasticity/impl/MajorSymmetric4Tensor.h>
#include <Elasticity/impl/Strain.h>
#include <sofa/core/objectmodel/BaseObject.h>

#if !defined(ELASTICITY_COMPONENT_HYPERELASTIC_MATERIAL_CPP)
#include <sofa/defaulttype/VecTypes.h>
#endif

namespace elasticity
{

template<class TDataTypes>
class HyperelasticMaterial : public virtual sofa::core::objectmodel::BaseObject
{
public:
    SOFA_ABSTRACT_CLASS(HyperelasticMaterial<TDataTypes>, sofa::core::objectmodel::BaseObject);
    using DataTypes = TDataTypes;

protected:
    using Real = sofa::Real_t<DataTypes>;

    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;

    using DeformationGradient = sofa::type::Mat<spatial_dimensions, spatial_dimensions, Real>;
    using RightCauchyGreenTensor = sofa::type::Mat<spatial_dimensions, spatial_dimensions, Real>;
    using StressTensor = sofa::type::Mat<spatial_dimensions, spatial_dimensions, Real>;
    using ElasticityTensor = elasticity::FullySymmetric4Tensor<DataTypes>;
    using TangentModulus = elasticity::MajorSymmetric4Tensor<DataTypes>;

public:
    void init() override;

    /**
     * Computes the First Piola-Kirchhoff stress tensor for a given deformation gradient.
     *
     * It corresponds to the derivative of the strain energy density function w.r.t. deformation
     * gradient.
     */
    virtual StressTensor firstPiolaKirchhoffStress(Strain<DataTypes>& strain) = 0;

    /**
     * Compute the jacobian of the first Piola-Kirchhoff stress tensor with respect to the
     * deformation gradient.
     *
     * It is called the material tangent modulus.
     */
    virtual TangentModulus materialTangentModulus(Strain<DataTypes>& strain) = 0;

};

#if !defined(ELASTICITY_COMPONENT_HYPERELASTIC_MATERIAL_CPP)
extern template class ELASTICITY_API HyperelasticMaterial<sofa::defaulttype::Vec1Types>;
extern template class ELASTICITY_API HyperelasticMaterial<sofa::defaulttype::Vec2Types>;
extern template class ELASTICITY_API HyperelasticMaterial<sofa::defaulttype::Vec3Types>;
#endif

}

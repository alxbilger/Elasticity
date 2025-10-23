#pragma once

#include <Elasticity/config.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <Elasticity/impl/Tensor.h>

#if !defined(ELASTICITY_COMPONENT_HYPERELASTIC_MATERIAL_CPP)
#include <sofa/defaulttype/VecTypes.h>
#endif

namespace elasticity
{

template<class DataTypes>
class HyperelasticMaterial : public virtual sofa::core::objectmodel::BaseObject
{
public:
    SOFA_ABSTRACT_CLASS(HyperelasticMaterial, sofa::core::objectmodel::BaseObject);

protected:
    using Real = sofa::Real_t<DataTypes>;

    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;

    using DeformationGradient = sofa::type::Mat<spatial_dimensions, spatial_dimensions, Real>;
    using StressTensor = sofa::type::Mat<spatial_dimensions, spatial_dimensions, Real>;
    using StressJacobian = elasticity::Tensor<Real,
        spatial_dimensions, spatial_dimensions, spatial_dimensions, spatial_dimensions>;

    static Real kroneckerDelta(std::size_t i, std::size_t j)
    {
        return static_cast<Real>(i == j);
    }

public:
    void init() override;

    /**
     * Computes the First Piola-Kirchhoff stress tensor for a given deformation gradient.
     *
     * It corresponds to the derivative of the strain energy density function w.r.t. deformation
     * gradient.
     */
    virtual StressTensor firstPiolaKirchhoffStress(const DeformationGradient& F);

    /**
     * Compute the jacobian of the first Piola-Kirchhoff stress tensor with respect to the
     * deformation gradient.
     *
     * It is called the material tangent modulus.
     */
    virtual StressJacobian materialTangentModulus(const DeformationGradient& F);

protected:

    virtual StressTensor secondPiolaKirchhoffStress(const DeformationGradient& F) { return {}; }
    virtual StressJacobian elasticityTensor(const DeformationGradient& F) { return {}; }

    /**
     * Compute the first Cauchy-Green invariant from the deformation gradient
     */
    static Real invariant1(const DeformationGradient& F);

    /**
     * Compute the second Cauchy-Green invariant from the deformation gradient
     */
    static Real invariant2(const DeformationGradient& F);

    /**
     * Compute the third Cauchy-Green invariant from the deformation gradient
     */
    static Real invariant3(const DeformationGradient& F);

};

#if !defined(ELASTICITY_COMPONENT_HYPERELASTIC_MATERIAL_CPP)
extern template class ELASTICITY_API HyperelasticMaterial<sofa::defaulttype::Vec1Types>;
extern template class ELASTICITY_API HyperelasticMaterial<sofa::defaulttype::Vec2Types>;
extern template class ELASTICITY_API HyperelasticMaterial<sofa::defaulttype::Vec3Types>;
#endif

}

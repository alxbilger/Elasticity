#pragma once

#include <Elasticity/config.h>
#include <sofa/core/objectmodel/BaseObject.h>

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
    using StressJacobian = sofa::type::Mat<
        spatial_dimensions * spatial_dimensions,
        spatial_dimensions * spatial_dimensions, Real>;

public:
    void init() override;

    /**
     * Computes the First Piola-Kirchhoff stress tensor for a given deformation gradient.
     */
    virtual StressTensor firstPiolaKirchhoffStress(const DeformationGradient& F) = 0;

    /**
     * Compute the jacobian of the first Piola-Kirchhoff stress tensor with respect to the deformation
     * gradient. The resulting 4th-order tensor must be flattened as a d^2 x d^2 matrix.
     */
    virtual StressJacobian jacobianFirstPiolaKirchhoffStress() = 0;
};

#if !defined(ELASTICITY_COMPONENT_HYPERELASTIC_MATERIAL_CPP)
extern template class ELASTICITY_API HyperelasticMaterial<sofa::defaulttype::Vec1Types>;
extern template class ELASTICITY_API HyperelasticMaterial<sofa::defaulttype::Vec2Types>;
extern template class ELASTICITY_API HyperelasticMaterial<sofa::defaulttype::Vec3Types>;
#endif

}

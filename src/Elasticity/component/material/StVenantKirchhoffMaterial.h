#pragma once

#include <Elasticity/config.h>
#include <Elasticity/component/HyperelasticMaterial.h>
#include <Elasticity/component/LinearMechanicalParametersComponent.h>

#if !defined(ELASTICITY_COMPONENT_MATERIAL_STVENANTKIRCHHOFFMATERIAL_CPP)
#include <sofa/defaulttype/VecTypes.h>
#endif

namespace elasticity
{

/**
 * @class StVenantKirchhoffMaterial
 * @brief Represents the St. Venant-Kirchhoff material model for hyperelastic materials.
 *
 * The St. Venant-Kirchhoff material is a simple model used to describe the stress-strain behavior
 * of isotropic hyperelastic materials. It is based on the linear elastic theory but is extended
 * to handle large deformations using the Green-Lagrange strain tensor and the Second
 * Piola-Kirchhoff stress tensor.
 *
 * This material model is only valid for cases where the deformation is relatively small, despite
 * being nonlinear in nature. Its application is mainly limited due to its inability to correctly
 * predict behavior under significant strain, as it does not accurately represent material
 * nonlinearity under large deformation.
 *
 * The material model is defined using two parameters:
 * - The Young's modulus: Describes the material's stiffness.
 * - The Poisson's ratio: Represents the material's ability to undergo deformation in directions
 *   orthogonal to the applied stress.
 */
template <class DataTypes>
class StVenantKirchhoffMaterial:
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

    using HyperelasticMaterial<DataTypes>::kroneckerDelta;

protected:
    StressTensor secondPiolaKirchhoffStress(const DeformationGradient& F) override;

    StressJacobian elasticityTensor(const DeformationGradient& F) override;
};

#if !defined(ELASTICITY_COMPONENT_MATERIAL_STVENANTKIRCHHOFFMATERIAL_CPP)
extern template class ELASTICITY_API StVenantKirchhoffMaterial<sofa::defaulttype::Vec1Types>;
extern template class ELASTICITY_API StVenantKirchhoffMaterial<sofa::defaulttype::Vec2Types>;
extern template class ELASTICITY_API StVenantKirchhoffMaterial<sofa::defaulttype::Vec3Types>;
#endif
}

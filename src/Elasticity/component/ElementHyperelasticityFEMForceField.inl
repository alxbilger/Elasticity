#pragma once

#include <Elasticity/component/ElementHyperelasticityFEMForceField.h>
#include <Elasticity/impl/ElasticityTensor.h>
#include <Elasticity/impl/MatrixTools.h>
#include <Elasticity/impl/VectorTools.h>
#include <sofa/core/ObjectFactory.h>

namespace elasticity
{

template <class DataTypes, class ElementType>
void ElementHyperelasticityFEMForceField<DataTypes, ElementType>::init()
{
    sofa::core::behavior::ForceField<DataTypes>::init();

    if (!this->isComponentStateInvalid())
    {
        this->validateTopology();
    }

    if (!this->isComponentStateInvalid())
    {
        this->validateMaterial();
    }

    if (!this->isComponentStateInvalid())
    {
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Valid);
    }
}

template <class DataTypes, class ElementType>
void ElementHyperelasticityFEMForceField<DataTypes, ElementType>::addForce(
    const sofa::core::MechanicalParams* mparams, DataVecDeriv& f, const DataVecCoord& x,
    const DataVecDeriv& v)
{
    SOFA_UNUSED(mparams);
    SOFA_UNUSED(v);

    auto forceAccessor = sofa::helper::getWriteOnlyAccessor(f);
    auto positionAccessor = sofa::helper::getReadAccessor(x);
    auto restPositionAccessor = this->mstate->readRestPositions();

    if (l_topology == nullptr) return;
    if (l_material == nullptr) return;

    const auto& elements = FiniteElement::getElementSequence(*l_topology);

    for (const auto& element : elements)
    {
        const std::array<Coord, NumberOfNodesInElement> elementNodesCoordinates = extractNodesVectorFromGlobalVector(element, positionAccessor.ref());
        const std::array<Coord, NumberOfNodesInElement> elementNodesRestCoordinates = extractNodesVectorFromGlobalVector(element, restPositionAccessor.ref());

        for (const auto& [quadraturePoint, weight] : FiniteElement::quadraturePoints())
        {
            // gradient of shape functions in the reference element evaluated at the quadrature point
            const sofa::type::Mat<NumberOfNodesInElement, ElementDimension, Real> dN_dq_ref =
                FiniteElement::gradientShapeFunctions(quadraturePoint);

            // jacobian of the mapping from the reference space to the physical space, evaluated at the
            // quadrature point
            sofa::type::Mat<spatial_dimensions, ElementDimension, Real> J_q;
            for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
                J_q += sofa::type::dyad(elementNodesCoordinates[i], dN_dq_ref[i]);

            // jacobian of the mapping from the reference space to the rest physical space, evaluated at the
            // quadrature point
            sofa::type::Mat<spatial_dimensions, ElementDimension, Real> J_Q;
            for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
                J_Q += sofa::type::dyad(elementNodesRestCoordinates[i], dN_dq_ref[i]);

            const auto detJ_Q = elasticity::determinant(J_Q);

            const sofa::type::Mat<ElementDimension, spatial_dimensions, Real> J_Q_inv = elasticity::inverse(J_Q);

            // gradient of the shape functions in the physical element evaluated at the quadrature point
            sofa::type::Mat<NumberOfNodesInElement, spatial_dimensions, Real> dN_dQ(sofa::type::NOINIT);
            for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
                dN_dQ[i] = J_Q_inv.transposed() * dN_dq_ref[i];

            // both ways to compute the deformation gradient are equivalent
            const DeformationGradient F = computeDeformationGradient(J_q, J_Q_inv);
            // const DeformationGradient F = computeDeformationGradient2(elementNodesCoordinates, dN_dQ);

            const auto detF = elasticity::determinant(F);
            if (detF < 0)
            {
                const auto detJ_q = elasticity::determinant(J_q);
                msg_error("FEM") << "Element inversion detected (detF = " << detF << " < 0, " <<
                    " detJ_q = " << detJ_q << ", detJ_Q = " << elasticity::determinant(J_Q) << ")";
            }

            const auto P = l_material->firstPiolaKirchhoffStress(F);

            for (sofa::Index i = 0; i < NumberOfNodesInElement; ++i)
            {
                forceAccessor[element[i]] += (-detJ_Q * weight) * P * dN_dQ[i];
            }
        }
    }

    // invalidate the Hessian, so it will be computed the next time it is necessary
    m_isHessianValid = false;
}

template <class DataTypes, class ElementType>
void ElementHyperelasticityFEMForceField<DataTypes, ElementType>::addDForce(
    const sofa::core::MechanicalParams* mparams, DataVecDeriv& df, const DataVecDeriv& dx)
{
    if (l_topology == nullptr) return;
    if (l_material == nullptr) return;

    auto dfAccessor = sofa::helper::getWriteAccessor(df);
    auto dxAccessor = sofa::helper::getReadAccessor(dx);
    dfAccessor.resize(dxAccessor.size());

    if (!m_isHessianValid)
    {
        computeHessian();
    }

    const auto& elements = FiniteElement::getElementSequence(*l_topology);

    const Real kFactor = static_cast<Real>(sofa::core::mechanicalparams::kFactorIncludingRayleighDamping(
            mparams, this->rayleighStiffness.getValue()));

    for (const auto& element : elements)
    {
    }
}

template <class DataTypes, class ElementType>
void ElementHyperelasticityFEMForceField<DataTypes, ElementType>::buildStiffnessMatrix(
    sofa::core::behavior::StiffnessMatrix* matrix)
{
    if (!m_isHessianValid)
    {
        computeHessian();
    }
}

template <class DataTypes, class ElementType>
SReal ElementHyperelasticityFEMForceField<DataTypes, ElementType>::getPotentialEnergy(
    const sofa::core::MechanicalParams*, const DataVecCoord& x) const
{
    return 0;
}

template <class DataTypes, class ElementType>
void ElementHyperelasticityFEMForceField<DataTypes, ElementType>::validateMaterial()
{
    if (l_material.empty())
    {
        msg_info() << "Link to a valid material should be set to ensure right behavior. First "
                      "material found in current context will be used.";
        l_material.set(this->getContext()->get<HyperelasticMaterial<DataTypes>>());
    }

    if (l_material == nullptr)
    {
        msg_error() << "No material component found at path: '" << this->l_material.getLinkedPath()
                    << "', nor in current context: " << this->getContext()->name
                    << ". Object must have a material. "
                    << "The list of available material components is: "
                    << sofa::core::ObjectFactory::getInstance()
                           ->listClassesDerivedFrom<HyperelasticMaterial<DataTypes>>();
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
    }
}

template <class DataTypes, class ElementType>
void ElementHyperelasticityFEMForceField<DataTypes, ElementType>::computeHessian()
{
    m_isHessianValid = true;
}

template <class DataTypes, class ElementType>
auto ElementHyperelasticityFEMForceField<DataTypes, ElementType>::computeDeformationGradient(
    const sofa::type::Mat<spatial_dimensions, ElementDimension, Real>& J_q,
    const sofa::type::Mat<ElementDimension, spatial_dimensions, Real>& J_Q_inv) -> DeformationGradient
{
    return J_q * J_Q_inv;
}

template <class DataTypes, class ElementType>
auto ElementHyperelasticityFEMForceField<DataTypes, ElementType>::computeDeformationGradient2(
    const std::array<Coord, NumberOfNodesInElement>& elementNodesCoordinates,
    const sofa::type::Mat<NumberOfNodesInElement, spatial_dimensions, Real>& dN_dQ)  -> DeformationGradient
{
    DeformationGradient F;

    for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
        F += sofa::type::dyad(elementNodesCoordinates[i], dN_dQ[i]);

    return F;
}

}  // namespace elasticity

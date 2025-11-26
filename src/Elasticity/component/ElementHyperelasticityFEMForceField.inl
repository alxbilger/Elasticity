#pragma once

#include <Elasticity/component/ElementHyperelasticityFEMForceField.h>
#include <Elasticity/impl/MatrixTools.h>
#include <Elasticity/impl/Strain.h>
#include <Elasticity/impl/VecView.h>
#include <Elasticity/impl/VectorTools.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/core/behavior/BaseLocalForceFieldMatrix.h>
#include <sofa/helper/ScopedAdvancedTimer.h>

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

    if (this->isComponentStateInvalid())
        return;

    auto forceAccessor = sofa::helper::getWriteOnlyAccessor(f);
    auto positionAccessor = sofa::helper::getReadAccessor(x);

    if (l_topology == nullptr) return;
    if (l_material == nullptr) return;

    m_coordinates = &positionAccessor.ref();

    const auto& elements = FiniteElement::getElementSequence(*l_topology);

    if (m_precomputedData.size() != elements.size())
    {
        precomputeData();
    }

    std::size_t elementIndex = 0;
    for (const auto& element : elements)
    {
        const std::array<Coord, NumberOfNodesInElement> elementNodesCoordinates = extractNodesVectorFromGlobalVector(element, positionAccessor.ref());

        static constexpr auto quadraturePoints = FiniteElement::quadraturePoints();
        static constexpr auto gradients = FiniteElementHelper<ElementType, DataTypes>::gradientShapeFunctionAtQuadraturePoints();

        for (sofa::Size q = 0; q < NumberOfQuadraturePoints; ++q)
        {
            const auto& weight = quadraturePoints[q].second;
            const PrecomputedData& precomputedData = m_precomputedData[elementIndex][q];

            // gradient of shape functions in the reference element evaluated at the quadrature point
            const sofa::type::Mat<NumberOfNodesInElement, TopologicalDimension, Real>& dN_dq_ref = gradients[q];

            // jacobian of the mapping from the reference space to the physical space, evaluated at the
            // quadrature point
            const auto J_q = FiniteElementHelper<ElementType, DataTypes>::jacobianFromReferenceToPhysical(elementNodesCoordinates, dN_dq_ref);

            const auto detJ_Q = precomputedData.detJacobian;

            const sofa::type::Mat<TopologicalDimension, spatial_dimensions, Real>& J_Q_inv = precomputedData.jacobianInv;

            // gradient of the shape functions in the physical element evaluated at the quadrature point
            const sofa::type::Mat<NumberOfNodesInElement, spatial_dimensions, Real>& dN_dQ = precomputedData.dN_dQ;

            // both ways to compute the deformation gradient are equivalent
            const DeformationGradient F = computeDeformationGradient(J_q, J_Q_inv);
            // const DeformationGradient F = computeDeformationGradient2(elementNodesCoordinates, dN_dQ);

            Strain<DataTypes> strain(deformationGradient, F);

            const auto detF = strain.getDeterminantDeformationGradient();
            if (detF < 0)
            {
                const auto detJ_q = elasticity::determinant(J_q);
                msg_error("FEM") << "Element inversion detected (detF = " << detF << " < 0, " <<
                    " detJ_q = " << detJ_q << ", detJ_Q = " << detJ_Q << ")";
            }

            const auto P = l_material->firstPiolaKirchhoffStress(strain);

            for (sofa::Index i = 0; i < NumberOfNodesInElement; ++i)
            {
                forceAccessor[element[i]] += (-detJ_Q * weight) * P * dN_dQ[i];
            }
        }
        ++elementIndex;
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

    if (this->isComponentStateInvalid())
        return;

    auto dfAccessor = sofa::helper::getWriteAccessor(df);
    auto dxAccessor = sofa::helper::getReadAccessor(dx);
    dfAccessor.resize(dxAccessor.size());

    if (!m_isHessianValid)
    {
        computeHessian(*m_coordinates);
    }

    const auto& elements = FiniteElement::getElementSequence(*l_topology);

    const Real kFactor = static_cast<Real>(sofa::core::mechanicalparams::kFactorIncludingRayleighDamping(
            mparams, this->rayleighStiffness.getValue()));

    auto elementStiffnessIt = m_elementStiffness.begin();
    for (const auto& element : elements)
    {
        sofa::type::Vec<NumberOfDofsInElement, Real> element_dx(sofa::type::NOINIT);

        for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
        {
            VecView<spatial_dimensions, Real> node_dx(element_dx, i * spatial_dimensions);
            node_dx = dxAccessor[element[i]];
        }

        const auto& stiffnessMatrix = *elementStiffnessIt++;
        auto dForce = (kFactor) * (stiffnessMatrix * element_dx);

        for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
        {
            VecView<spatial_dimensions, Real> nodedForce(dForce, i * spatial_dimensions);
            dfAccessor[element[i]] += nodedForce.toVec();
        }
    }
}

template <class DataTypes, class ElementType>
void ElementHyperelasticityFEMForceField<DataTypes, ElementType>::buildStiffnessMatrix(
    sofa::core::behavior::StiffnessMatrix* matrix)
{
    if (this->isComponentStateInvalid())
        return;

    auto dfdx = matrix->getForceDerivativeIn(this->mstate)
        .withRespectToPositionsIn(this->mstate);

    if (!m_isHessianValid)
    {
        computeHessian(*m_coordinates);
    }

    sofa::type::Mat<spatial_dimensions, spatial_dimensions, Real> localMatrix(sofa::type::NOINIT);

    const auto& elements = FiniteElement::getElementSequence(*l_topology);
    auto elementStiffnessIt = m_elementStiffness.begin();
    for (const auto& element : elements)
    {
        const auto& stiffnessMatrix = *elementStiffnessIt++;

        for (sofa::Index n1 = 0; n1 < NumberOfNodesInElement; ++n1)
        {
            for (sofa::Index n2 = 0; n2 < NumberOfNodesInElement; ++n2)
            {
                stiffnessMatrix.getsub(spatial_dimensions * n1, spatial_dimensions * n2, localMatrix); //extract the submatrix corresponding to the coupling of nodes n1 and n2
                dfdx(element[n1] * spatial_dimensions, element[n2] * spatial_dimensions) += localMatrix;
            }
        }
    }
}

template <class DataTypes, class ElementType>
SReal ElementHyperelasticityFEMForceField<DataTypes, ElementType>::getPotentialEnergy(
    const sofa::core::MechanicalParams*, const DataVecCoord& x) const
{
    return 0;
}

template <class TDataTypes, class TElementType>
void ElementHyperelasticityFEMForceField<TDataTypes, TElementType>::addKToMatrix(
    sofa::linearalgebra::BaseMatrix* matrix, SReal kFact, unsigned& offset)
{
    if (this->isComponentStateInvalid())
        return;

    if (!m_isHessianValid)
    {
        computeHessian(*m_coordinates);
    }

    using LocalMatType = sofa::type::Mat<spatial_dimensions, spatial_dimensions, Real>;
    LocalMatType localMatrix{sofa::type::NOINIT};

    const auto& elements = FiniteElement::getElementSequence(*l_topology);
    auto elementStiffnessIt = m_elementStiffness.begin();
    for (const auto& element : elements)
    {
        const auto& stiffnessMatrix = *elementStiffnessIt++;

        for (sofa::Index n1 = 0; n1 < NumberOfNodesInElement; ++n1)
        {
            for (sofa::Index n2 = 0; n2 < NumberOfNodesInElement; ++n2)
            {
                stiffnessMatrix.getsub(spatial_dimensions * n1, spatial_dimensions * n2, localMatrix); //extract the submatrix corresponding to the coupling of nodes n1 and n2
                const auto value = (-static_cast<Real>(kFact)) * static_cast<ScalarOrMatrix<LocalMatType>>(localMatrix);
                matrix->add(
                   offset + element[n1] * spatial_dimensions,
                   offset + element[n2] * spatial_dimensions, value);
            }
        }
    }
}

template <class DataTypes, class ElementType>
void ElementHyperelasticityFEMForceField<DataTypes, ElementType>::validateMaterial()
{
    if (l_material.empty())
    {
        msg_info() << "Link to a valid material should be set to ensure right behavior. First "
                      "material found in current context will be used.";
        l_material.set(this->getContext()->template get<HyperelasticMaterial<DataTypes>>());
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
void ElementHyperelasticityFEMForceField<DataTypes, ElementType>::computeHessian(const VecCoord& coordinates)
{
    if (l_topology == nullptr) return;
    if (l_material == nullptr) return;

    SCOPED_TIMER("ComputeHessian");

    const auto& elements = FiniteElement::getElementSequence(*l_topology);

    if (m_precomputedData.size() != elements.size())
    {
        precomputeData();
    }

    m_elementStiffness.clear();
    m_elementStiffness.resize(elements.size());

    auto elementStiffnessIt = m_elementStiffness.begin();

    std::size_t elementIndex = 0;
    for (const auto& element : elements)
    {
        SCOPED_TIMER_TR("Element");
        const std::array<Coord, NumberOfNodesInElement> elementNodesCoordinates = extractNodesVectorFromGlobalVector(element, this->mstate->readPositions().ref());

        ElementStiffness& K = *elementStiffnessIt++;

        static constexpr auto quadraturePoints = FiniteElement::quadraturePoints();
        static constexpr auto gradients = FiniteElementHelper<ElementType, DataTypes>::gradientShapeFunctionAtQuadraturePoints();

        for (sofa::Size q = 0; q < NumberOfQuadraturePoints; ++q)
        {
            const auto& weight = quadraturePoints[q].second;
            const PrecomputedData& precomputedData = m_precomputedData[elementIndex][q];

            // gradient of shape functions in the reference element evaluated at the quadrature point
            const sofa::type::Mat<NumberOfNodesInElement, TopologicalDimension, Real>& dN_dq_ref = gradients[q];

            // jacobian of the mapping from the reference space to the physical space, evaluated at the
            // quadrature point
            const auto J_q = FiniteElementHelper<ElementType, DataTypes>::jacobianFromReferenceToPhysical(elementNodesCoordinates, dN_dq_ref);

            const auto detJ_Q = precomputedData.detJacobian;

            const sofa::type::Mat<TopologicalDimension, spatial_dimensions, Real>& J_Q_inv = precomputedData.jacobianInv;

            // gradient of the shape functions in the physical element evaluated at the quadrature point
            const sofa::type::Mat<NumberOfNodesInElement, spatial_dimensions, Real>& dN_dQ = precomputedData.dN_dQ;

            // both ways to compute the deformation gradient are equivalent
            const DeformationGradient F = computeDeformationGradient(J_q, J_Q_inv);
            // const DeformationGradient F = computeDeformationGradient2(elementNodesCoordinates, dN_dQ);

            Strain<DataTypes> strain(deformationGradient, F);

            // derivative of first Piola-Kirchhoff stress tensor with respect to deformation gradient
            const auto dPdF = [&]()
            {
                SCOPED_TIMER_VARNAME_TR(dPdFTimer, "dPdF");
                return l_material->materialTangentModulus(strain);
            }();

            const auto factor = -detJ_Q * weight;

            SCOPED_TIMER_VARNAME_TR(tensorTimer, "Tensor");
            for (sofa::Size element_i = 0; element_i < NumberOfNodesInElement; ++element_i)
            {
                for (sofa::Size element_j = 0; element_j < NumberOfNodesInElement; ++element_j)
                {
                    for (sofa::Size dimension_i = 0; dimension_i < spatial_dimensions; ++dimension_i)
                    {
                        for (sofa::Size dimension_j = 0; dimension_j < spatial_dimensions; ++dimension_j)
                        {
                            auto& k = K(element_i * spatial_dimensions + dimension_i, element_j * spatial_dimensions + dimension_j);
                            for (sofa::Size i = 0; i < spatial_dimensions; ++i)
                            {
                                for (sofa::Size j = 0; j < spatial_dimensions; ++j)
                                {
                                     k += factor * dN_dQ[element_i][i] * dPdF(dimension_i, i, dimension_j, j) * dN_dQ[element_j][j];
                                }
                            }
                        }
                    }
                }
            }
        }
        ++elementIndex;
    }

    m_isHessianValid = true;
}

template <class DataTypes, class ElementType>
auto ElementHyperelasticityFEMForceField<DataTypes, ElementType>::computeDeformationGradient(
    const sofa::type::Mat<spatial_dimensions, TopologicalDimension, Real>& J_q,
    const sofa::type::Mat<TopologicalDimension, spatial_dimensions, Real>& J_Q_inv) -> DeformationGradient
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

template <class TDataTypes, class TElementType>
void ElementHyperelasticityFEMForceField<TDataTypes, TElementType>::precomputeData()
{
    if (l_topology == nullptr) return;

    auto restPositionAccessor = this->mstate->readRestPositions();
    const auto& restPosition = restPositionAccessor.ref();

    const auto& elements = FiniteElement::getElementSequence(*l_topology);
    m_precomputedData.resize(elements.size());

    static constexpr auto gradients = FiniteElementHelper<TElementType, TDataTypes>::gradientShapeFunctionAtQuadraturePoints();

    for (std::size_t i = 0; i < elements.size(); ++i)
    {
        const auto& element = elements[i];
        const std::array<Coord, NumberOfNodesInElement> elementNodesRestCoordinates = extractNodesVectorFromGlobalVector(element, restPosition);

        for (std::size_t j = 0; j < NumberOfQuadraturePoints; ++j)
        {
            // gradient of shape functions in the reference element evaluated at the quadrature point
            const sofa::type::Mat<NumberOfNodesInElement, TopologicalDimension, Real>& dN_dq_ref = gradients[j];

            PrecomputedData& data = m_precomputedData[i][j];
            data.jacobian = FiniteElementHelper<TElementType, TDataTypes>::jacobianFromReferenceToPhysical(elementNodesRestCoordinates, dN_dq_ref);
            data.jacobianInv = elasticity::inverse(data.jacobian);
            data.detJacobian = elasticity::determinant(data.jacobian);

            for (sofa::Size n = 0; n < NumberOfNodesInElement; ++n)
            {
                data.dN_dQ[n] = data.jacobianInv.multTranspose(dN_dq_ref[n]);
            }
        }
    }
}

}  // namespace elasticity

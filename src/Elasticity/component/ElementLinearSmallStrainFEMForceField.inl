#pragma once
#include <Elasticity/component/ElementLinearSmallStrainFEMForceField.h>
#include <Elasticity/impl/VecView.h>
#include <Elasticity/impl/VectorTools.h>

#include <sofa/core/behavior/ForceField.inl>

namespace elasticity
{

template <class DataTypes, class ElementType>
void ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::init()
{
    sofa::core::behavior::ForceField<DataTypes>::init();

    if (!this->isComponentStateInvalid())
    {
        this->validateTopology();
    }

    if (!this->isComponentStateInvalid())
    {
        precomputeElementStiffness();
    }

    if (!this->isComponentStateInvalid())
    {
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Valid);
    }
}

template <class DataTypes, class ElementType>
void ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::addForce(
    const sofa::core::MechanicalParams* mparams, DataVecDeriv& f, const DataVecCoord& x,
    const DataVecDeriv& v)
{
    SOFA_UNUSED(mparams);
    SOFA_UNUSED(v);

    auto forceAccessor = sofa::helper::getWriteOnlyAccessor(f);
    auto positionAccessor = sofa::helper::getReadAccessor(x);
    auto restPositionAccessor = this->mstate->readRestPositions();

    if (l_topology == nullptr) return;

    const auto& elements = FiniteElement::getElementSequence(*l_topology);

    for (sofa::Size i = 0; i < elements.size(); ++i)
    {
        const auto& element = elements[i];
        const ElementStiffness& stiffnessMatrix = m_elementStiffness[i];

        const std::array<Coord, NumberOfNodesInElement> elementNodesCoordinates = extractNodesVectorFromGlobalVector(element, positionAccessor.ref());
        const std::array<Coord, NumberOfNodesInElement> restElementNodesCoordinates = extractNodesVectorFromGlobalVector(element, restPositionAccessor.ref());

        const ElementDisplacement displacement = computeElementDisplacement(elementNodesCoordinates, restElementNodesCoordinates);

        sofa::type::Vec<NumberOfDofsInElement, Real> elementForce = stiffnessMatrix * displacement;

        for (sofa::Size j = 0; j < NumberOfNodesInElement; ++j)
        {
            VecView<spatial_dimensions, Real> nodeForce(elementForce, j * spatial_dimensions);
            forceAccessor[element[j]] += -nodeForce;
        }
    }
}

template <class DataTypes, class ElementType>
void ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::addDForce(
    const sofa::core::MechanicalParams* mparams, DataVecDeriv& df, const DataVecDeriv& dx)
{
    auto dfAccessor = sofa::helper::getWriteAccessor(df);
    auto dxAccessor = sofa::helper::getReadAccessor(dx);
    dfAccessor.resize(dxAccessor.size());

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
        auto dForce = (-kFactor) * (stiffnessMatrix * element_dx);

        for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
        {
            VecView<spatial_dimensions, Real> nodedForce(dForce, i * spatial_dimensions);
            dfAccessor[element[i]] += nodedForce.toVec();
        }
    }
}

template <class DataTypes, class ElementType>
void ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::buildStiffnessMatrix(
    sofa::core::behavior::StiffnessMatrix* matrix)
{
    auto dfdx = matrix->getForceDerivativeIn(this->mstate)
        .withRespectToPositionsIn(this->mstate);

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
                dfdx(element[n1] * spatial_dimensions, element[n2] * spatial_dimensions) += -localMatrix;
            }
        }
    }
}

template <class DataTypes, class ElementType>
SReal ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::getPotentialEnergy(
    const sofa::core::MechanicalParams*, const DataVecCoord& x) const
{
    return 0;
}

template <class DataTypes, class ElementType>
void ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::addKToMatrix(
    sofa::linearalgebra::BaseMatrix* matrix, SReal kFact, unsigned& offset)
{
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
                matrix->add(
                    offset + element[n1] * spatial_dimensions,
                    offset + element[n2] * spatial_dimensions, -kFact * localMatrix);
            }
        }
    }
}

template <class DataTypes, class ElementType>
void ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::precomputeElementStiffness()
{
    const auto youngModulus = this->d_youngModulus.getValue();
    const auto poissonRatio = this->d_poissonRatio.getValue();

    auto restPositionAccessor = this->mstate->readRestPositions();

    m_elasticityTensor = makeIsotropicElasticityTensor<DataTypes>(youngModulus, poissonRatio);

    m_elementStiffness.clear();

    const auto& elements = FiniteElement::getElementSequence(*l_topology);
    m_elementStiffness.reserve(elements.size());

    for (const auto& element : elements)
    {
        const std::array<Coord, NumberOfNodesInElement> nodesCoordinates = extractNodesVectorFromGlobalVector(element, restPositionAccessor.ref());
        ElementStiffness K = integrate<DataTypes, ElementType>(nodesCoordinates, m_elasticityTensor);
        m_elementStiffness.push_back(K);
    }

    // precompute strain-displacement tensors at the nodes positions
    m_strainDisplacement.clear();
    m_strainDisplacement.resize(elements.size());

    static const std::array<ReferenceCoord, NumberOfNodesInElement>& referenceElementNodes =
        FiniteElement::referenceElementNodes;

    for (sofa::Size i = 0; i < elements.size(); ++i)
    {
        const auto& element = elements[i];
        const auto nodesCoordinates = extractNodesVectorFromGlobalVector(element, restPositionAccessor.ref());
        for (sofa::Size j = 0; j < NumberOfNodesInElement; ++j)
        {
            const ReferenceCoord& x = referenceElementNodes[j];

            // gradient of shape functions in the reference element evaluated at the quadrature point
            const sofa::type::Mat<NumberOfNodesInElement, ElementDimension, Real> dN_dq_ref =
                FiniteElement::gradientShapeFunctions(x);

            // jacobian of the mapping from the reference space to the physical space, evaluated at the
            // quadrature point
            sofa::type::Mat<spatial_dimensions, ElementDimension, Real> jacobian;
            for (sofa::Size n = 0; n < NumberOfNodesInElement; ++n)
                jacobian += sofa::type::dyad(nodesCoordinates[n], dN_dq_ref[n]);

            const sofa::type::Mat<ElementDimension, spatial_dimensions, Real> J_inv =
                elasticity::inverse(jacobian);

            // gradient of the shape functions in the physical element evaluated at the quadrature point
            sofa::type::Mat<NumberOfNodesInElement, spatial_dimensions, Real> dN_dq(sofa::type::NOINIT);
            for (sofa::Size n = 0; n < NumberOfNodesInElement; ++n)
                dN_dq[n] = J_inv.transposed() * dN_dq_ref[n];

            const auto B = makeStrainDisplacement<DataTypes, ElementType>(dN_dq);

            m_strainDisplacement[i][j] = B;
        }
    }
}

template <class DataTypes, class ElementType>
auto ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::computeElementDisplacement(
    const std::array<Coord, NumberOfNodesInElement>& elementNodesCoordinates,
    const std::array<Coord, NumberOfNodesInElement>& restElementNodesCoordinates) -> ElementDisplacement
{
    ElementDisplacement displacement(sofa::type::NOINIT);
    for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
    {
        for (sofa::Size j = 0; j < spatial_dimensions; ++j)
        {
            displacement[i * spatial_dimensions + j] = elementNodesCoordinates[i][j] - restElementNodesCoordinates[i][j];
        }
    }
    return displacement;
}

}  // namespace elasticity

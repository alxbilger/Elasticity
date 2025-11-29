#pragma once
#include <Elasticity/component/ElementCorotationalFEMForceField.h>
#include <Elasticity/component/ElementLinearSmallStrainFEMForceField.inl>

namespace elasticity
{

template <class DataTypes, class ElementType>
void ElementCorotationalFEMForceField<DataTypes, ElementType>::addForce(
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

    m_rotations.resize(elements.size(), RotationMatrix::Identity());

    auto elementStiffnessIt = this->m_elementStiffness.begin();
    auto rotationMatrixIt = m_rotations.begin();

    for (const auto element : elements)
    {
        const std::array<Coord, NumberOfNodesInElement> elementNodesCoordinates =
            extractNodesVectorFromGlobalVector(element, positionAccessor.ref());
        const std::array<Coord, NumberOfNodesInElement> restElementNodesCoordinates =
            extractNodesVectorFromGlobalVector(element, restPositionAccessor.ref());

        auto& elementRotation = *rotationMatrixIt++;
        computeElementRotation(elementNodesCoordinates, restElementNodesCoordinates, elementRotation);

        const auto t = translation(elementNodesCoordinates);
        const auto t0 = translation(restElementNodesCoordinates);

        ElementDisplacement displacement(sofa::type::NOINIT);
        for (sofa::Size j = 0; j < NumberOfNodesInElement; ++j)
        {
            VecView<spatial_dimensions, Real> transformedDisplacement(displacement, j * spatial_dimensions);
            transformedDisplacement = elementRotation.multTranspose(elementNodesCoordinates[j] - t) - (restElementNodesCoordinates[j] - t0);
        }

        const ElementStiffness& stiffnessMatrix = *elementStiffnessIt++;
        sofa::type::Vec<NumberOfDofsInElement, Real> elementForce = stiffnessMatrix * displacement;

        for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
        {
            VecView<spatial_dimensions, Real> nodeForce(elementForce, i * spatial_dimensions);
            forceAccessor[element[i]] += -elementRotation * nodeForce;
        }
    }
}

template <class DataTypes, class ElementType>
void ElementCorotationalFEMForceField<DataTypes, ElementType>::addDForce(
    const sofa::core::MechanicalParams* mparams, DataVecDeriv& df, const DataVecDeriv& dx)
{
    auto dfAccessor = sofa::helper::getWriteAccessor(df);
    auto dxAccessor = sofa::helper::getReadAccessor(dx);
    dfAccessor.resize(dxAccessor.size());

    const auto& elements = FiniteElement::getElementSequence(*l_topology);

    const Real kFactor = static_cast<Real>(sofa::core::mechanicalparams::kFactorIncludingRayleighDamping(
        mparams, this->rayleighStiffness.getValue()));

    sofa::type::Vec<NumberOfDofsInElement, Real> dForce;

    for (sofa::Size e = 0; e < elements.size(); ++e)
    {
        const auto& element = elements[e];
        const auto& elementRotation = m_rotations[e];

        sofa::type::Vec<NumberOfDofsInElement, Real> element_dx(sofa::type::NOINIT);

        for (sofa::Size n = 0; n < NumberOfNodesInElement; ++n)
        {
            VecView<spatial_dimensions, Real> rotated_dx(element_dx, n * spatial_dimensions);
            rotated_dx = elementRotation.multTranspose(dxAccessor[element[n]]);
        }

        const auto& stiffnessMatrix = this->m_elementStiffness[e];
        dForce = (-kFactor) * (stiffnessMatrix * element_dx);

        for (sofa::Size n = 0; n < NumberOfNodesInElement; ++n)
        {
            VecView<spatial_dimensions, Real> nodedForce(dForce, n * spatial_dimensions);
            dfAccessor[element[n]] += elementRotation * nodedForce;
        }
    }
}

template <class DataTypes, class ElementType>
void ElementCorotationalFEMForceField<DataTypes, ElementType>::buildStiffnessMatrix(
    sofa::core::behavior::StiffnessMatrix* matrix)
{
    auto dfdx = matrix->getForceDerivativeIn(this->mstate)
        .withRespectToPositionsIn(this->mstate);

    sofa::type::Mat<spatial_dimensions, spatial_dimensions, Real> localMatrix(sofa::type::NOINIT);

    const auto& elements = FiniteElement::getElementSequence(*l_topology);
    auto elementStiffnessIt = this->m_elementStiffness.begin();
    auto rotationMatrixIt = m_rotations.begin();
    for (const auto& element : elements)
    {
        const auto& elementRotation = *rotationMatrixIt++;
        const auto elementRotation_T = elementRotation.transposed();

        const auto& stiffnessMatrix = *elementStiffnessIt++;

        for (sofa::Index n1 = 0; n1 < NumberOfNodesInElement; ++n1)
        {
            for (sofa::Index n2 = 0; n2 < NumberOfNodesInElement; ++n2)
            {
                stiffnessMatrix.getsub(spatial_dimensions * n1, spatial_dimensions * n2, localMatrix);  // extract the submatrix corresponding to the coupling of nodes n1 and n2
                dfdx(element[n1] * spatial_dimensions, element[n2] * spatial_dimensions) += -elementRotation * localMatrix * elementRotation_T;
            }
        }
    }
}

template <class DataTypes, class ElementType>
SReal ElementCorotationalFEMForceField<DataTypes, ElementType>::getPotentialEnergy(
    const sofa::core::MechanicalParams*, const DataVecCoord& x) const
{
    return 0;
}

template <class DataTypes, class ElementType>
void ElementCorotationalFEMForceField<DataTypes, ElementType>::computeElementRotation(
    const std::array<Coord, NumberOfNodesInElement>& nodesPosition,
    const std::array<Coord, NumberOfNodesInElement>& nodesRestPosition,
    RotationMatrix& rotationMatrix)
{
    const auto t = translation(nodesPosition);
    const auto t0 = translation(nodesRestPosition);

    sofa::type::Mat<NumberOfNodesInElement, spatial_dimensions, Real> P(sofa::type::NOINIT);
    sofa::type::Mat<NumberOfNodesInElement, spatial_dimensions, Real> Q(sofa::type::NOINIT);

    for (sofa::Size j = 0; j < NumberOfNodesInElement; ++j)
    {
        P[j] = nodesPosition[j] - t;
        Q[j] = nodesRestPosition[j] - t0;
    }

    const auto H = P.multTranspose(Q);

    sofa::helper::Decompose<Real>::polarDecomposition_stable(H, rotationMatrix);
}

template <class DataTypes, class ElementType>
auto ElementCorotationalFEMForceField<DataTypes, ElementType>::translation(
    const std::array<Coord, NumberOfNodesInElement>& nodes) const -> Coord
{
    return computeCentroid(nodes);
}

template <class DataTypes, class ElementType>
auto ElementCorotationalFEMForceField<DataTypes, ElementType>::computeCentroid(
    const std::array<Coord, NumberOfNodesInElement>& nodes) -> Coord
{
    Coord centroid;
    for (const auto node : nodes)
    {
        centroid += node;
    }
    centroid /= static_cast<Real>(NumberOfNodesInElement);
    return centroid;
}

}  // namespace elasticity

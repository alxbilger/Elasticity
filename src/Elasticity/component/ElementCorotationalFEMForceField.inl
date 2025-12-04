#pragma once
#include <Elasticity/component/ElementCorotationalFEMForceField.h>
#include <Elasticity/component/ElementLinearSmallStrainFEMForceField.inl>

namespace elasticity
{

template <class DataTypes, class ElementType>
void ElementCorotationalFEMForceField<DataTypes, ElementType>::addForce(
    const sofa::core::MechanicalParams* mparams,
    sofa::DataVecDeriv_t<DataTypes>& f,
    const sofa::DataVecCoord_t<DataTypes>& x,
    const sofa::DataVecDeriv_t<DataTypes>& v)
{
    SOFA_UNUSED(mparams);
    SOFA_UNUSED(v);

    auto forceAccessor = sofa::helper::getWriteOnlyAccessor(f);
    auto positionAccessor = sofa::helper::getReadAccessor(x);
    auto restPositionAccessor = this->mstate->readRestPositions();

    if (l_topology == nullptr) return;

    const auto& elements = trait::FiniteElement::getElementSequence(*l_topology);

    m_rotations.resize(elements.size(), RotationMatrix::Identity());

    auto elementStiffnessIt = this->m_elementStiffness.begin();
    auto rotationMatrixIt = m_rotations.begin();

    for (const auto element : elements)
    {
        const std::array<sofa::Coord_t<DataTypes>, trait::NumberOfNodesInElement> elementNodesCoordinates =
            extractNodesVectorFromGlobalVector(element, positionAccessor.ref());
        const std::array<sofa::Coord_t<DataTypes>, trait::NumberOfNodesInElement> restElementNodesCoordinates =
            extractNodesVectorFromGlobalVector(element, restPositionAccessor.ref());

        auto& elementRotation = *rotationMatrixIt++;
        computeElementRotation(elementNodesCoordinates, restElementNodesCoordinates, elementRotation);

        const auto t = translation(elementNodesCoordinates);
        const auto t0 = translation(restElementNodesCoordinates);

        typename trait::ElementDisplacement displacement(sofa::type::NOINIT);
        for (sofa::Size j = 0; j < trait::NumberOfNodesInElement; ++j)
        {
            VecView<trait::spatial_dimensions, sofa::Real_t<DataTypes>> transformedDisplacement(displacement, j * trait::spatial_dimensions);
            transformedDisplacement = elementRotation.multTranspose(elementNodesCoordinates[j] - t) - (restElementNodesCoordinates[j] - t0);
        }

        const auto& stiffnessMatrix = *elementStiffnessIt++;
        sofa::type::Vec<trait::NumberOfDofsInElement, sofa::Real_t<DataTypes>> elementForce = stiffnessMatrix * displacement;

        for (sofa::Size i = 0; i < trait::NumberOfNodesInElement; ++i)
        {
            VecView<trait::spatial_dimensions, sofa::Real_t<DataTypes>> nodeForce(elementForce, i * trait::spatial_dimensions);
            forceAccessor[element[i]] += -elementRotation * nodeForce;
        }
    }
}

template <class DataTypes, class ElementType>
void ElementCorotationalFEMForceField<DataTypes, ElementType>::addDForce(
    const sofa::core::MechanicalParams* mparams,
    sofa::DataVecDeriv_t<DataTypes>& df, const
    sofa::DataVecDeriv_t<DataTypes>& dx)
{
    auto dfAccessor = sofa::helper::getWriteAccessor(df);
    auto dxAccessor = sofa::helper::getReadAccessor(dx);
    dfAccessor.resize(dxAccessor.size());

    const auto& elements = trait::FiniteElement::getElementSequence(*l_topology);

    const sofa::Real_t<DataTypes> kFactor = static_cast<sofa::Real_t<DataTypes>>(sofa::core::mechanicalparams::kFactorIncludingRayleighDamping(
        mparams, this->rayleighStiffness.getValue()));

    sofa::type::Vec<trait::NumberOfDofsInElement, sofa::Real_t<DataTypes>> dForce;

    for (sofa::Size e = 0; e < elements.size(); ++e)
    {
        const auto& element = elements[e];
        const auto& elementRotation = m_rotations[e];

        sofa::type::Vec<trait::NumberOfDofsInElement, sofa::Real_t<DataTypes>> element_dx(sofa::type::NOINIT);

        for (sofa::Size n = 0; n < trait::NumberOfNodesInElement; ++n)
        {
            VecView<trait::spatial_dimensions, sofa::Real_t<DataTypes>> rotated_dx(element_dx, n * trait::spatial_dimensions);
            rotated_dx = elementRotation.multTranspose(dxAccessor[element[n]]);
        }

        const auto& stiffnessMatrix = this->m_elementStiffness[e];
        dForce = (-kFactor) * (stiffnessMatrix * element_dx);

        for (sofa::Size n = 0; n < trait::NumberOfNodesInElement; ++n)
        {
            VecView<trait::spatial_dimensions, sofa::Real_t<DataTypes>> nodedForce(dForce, n * trait::spatial_dimensions);
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

    sofa::type::Mat<trait::spatial_dimensions, trait::spatial_dimensions, sofa::Real_t<DataTypes>> localMatrix(sofa::type::NOINIT);

    const auto& elements = trait::FiniteElement::getElementSequence(*l_topology);
    auto elementStiffnessIt = this->m_elementStiffness.begin();
    auto rotationMatrixIt = m_rotations.begin();
    for (const auto& element : elements)
    {
        const auto& elementRotation = *rotationMatrixIt++;
        const auto elementRotation_T = elementRotation.transposed();

        const auto& stiffnessMatrix = *elementStiffnessIt++;

        for (sofa::Index n1 = 0; n1 < trait::NumberOfNodesInElement; ++n1)
        {
            for (sofa::Index n2 = 0; n2 < trait::NumberOfNodesInElement; ++n2)
            {
                stiffnessMatrix.getAssembledMatrix().getsub(trait::spatial_dimensions * n1, trait::spatial_dimensions * n2, localMatrix);  // extract the submatrix corresponding to the coupling of nodes n1 and n2
                dfdx(element[n1] * trait::spatial_dimensions, element[n2] * trait::spatial_dimensions) += -elementRotation * localMatrix * elementRotation_T;
            }
        }
    }
}

template <class DataTypes, class ElementType>
SReal ElementCorotationalFEMForceField<DataTypes, ElementType>::getPotentialEnergy(
    const sofa::core::MechanicalParams*,
    const sofa::DataVecCoord_t<DataTypes>& x) const
{
    return 0;
}

template <class DataTypes, class ElementType>
void ElementCorotationalFEMForceField<DataTypes, ElementType>::computeElementRotation(
    const std::array<sofa::Coord_t<DataTypes>, trait::NumberOfNodesInElement>& nodesPosition,
    const std::array<sofa::Coord_t<DataTypes>, trait::NumberOfNodesInElement>& nodesRestPosition,
    RotationMatrix& rotationMatrix)
{
    const auto t = translation(nodesPosition);
    const auto t0 = translation(nodesRestPosition);

    sofa::type::Mat<trait::NumberOfNodesInElement, trait::spatial_dimensions, sofa::Real_t<DataTypes>> P(sofa::type::NOINIT);
    sofa::type::Mat<trait::NumberOfNodesInElement, trait::spatial_dimensions, sofa::Real_t<DataTypes>> Q(sofa::type::NOINIT);

    for (sofa::Size j = 0; j < trait::NumberOfNodesInElement; ++j)
    {
        P[j] = nodesPosition[j] - t;
        Q[j] = nodesRestPosition[j] - t0;
    }

    const auto H = P.multTranspose(Q);

    sofa::helper::Decompose<sofa::Real_t<DataTypes>>::polarDecomposition_stable(H, rotationMatrix);
}

template <class DataTypes, class ElementType>
auto ElementCorotationalFEMForceField<DataTypes, ElementType>::translation(
    const std::array<sofa::Coord_t<DataTypes>, trait::NumberOfNodesInElement>& nodes) const -> sofa::Coord_t<DataTypes>
{
    return computeCentroid(nodes);
}

template <class DataTypes, class ElementType>
auto ElementCorotationalFEMForceField<DataTypes, ElementType>::computeCentroid(
    const std::array<sofa::Coord_t<DataTypes>, trait::NumberOfNodesInElement>& nodes) -> sofa::Coord_t<DataTypes>
{
    sofa::Coord_t<DataTypes> centroid;
    for (const auto node : nodes)
    {
        centroid += node;
    }
    centroid /= static_cast<sofa::Real_t<DataTypes>>(trait::NumberOfNodesInElement);
    return centroid;
}

}  // namespace elasticity

#pragma once
#include <Elasticity/CorotationalFEM.h>
#include <Elasticity/MatrixTools.h>
#include <Elasticity/VectorTools.h>
#include <sofa/helper/decompose.h>

namespace elasticity
{

template <class DataTypes, class ElementType>
void CorotationalFEM<DataTypes, ElementType>::addForce(VecDeriv& force, const VecCoord& position,
                                                       const VecCoord& restPosition)
{
    if (m_topology == nullptr) return;

    const auto& elements = FiniteElement::getElementSequence(*m_topology);

    Deriv nodeForce(sofa::type::NOINIT);

    m_rotations.resize(elements.size(), RotationMatrix::Identity());

    auto elementStiffnessIt = this->stiffnessMatrices().begin();
    auto rotationMatrixIt = m_rotations.begin();

    for (const auto element : elements)
    {
        const std::array<Coord, NumberOfNodesInElement> elementNodesCoordinates =
            extractNodesVectorFromGlobalVector(element, position);
        const std::array<Coord, NumberOfNodesInElement> restElementNodesCoordinates =
            extractNodesVectorFromGlobalVector(element, restPosition);

        auto& elementRotation = *rotationMatrixIt++;
        elementRotation =
            computeElementRotation(elementNodesCoordinates, restElementNodesCoordinates);

        const auto t = translation(elementNodesCoordinates);
        std::array<Coord, NumberOfNodesInElement> rotatedDisplacement;
        for (sofa::Size j = 0; j < NumberOfNodesInElement; ++j)
        {
            rotatedDisplacement[j] =
                elementRotation.transposed() * (elementNodesCoordinates[j] - t) -
                restElementNodesCoordinates[j];
        }

        ElementDisplacement displacement(sofa::type::NOINIT);
        for (sofa::Size j = 0; j < NumberOfNodesInElement; ++j)
        {
            for (sofa::Size k = 0; k < spatial_dimensions; ++k)
            {
                displacement[j * spatial_dimensions + k] = rotatedDisplacement[j][k];
            }
        }

        const ElementStiffness& stiffnessMatrix = *elementStiffnessIt++;
        const sofa::type::Vec<NumberOfDofsInElement, Real> elementForce =
            stiffnessMatrix * displacement;

        for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
        {
            elementForce.getsub(i * spatial_dimensions, nodeForce);
            force[element[i]] += -elementRotation * nodeForce;
        }
    }
}

template <class DataTypes, class ElementType>
void CorotationalFEM<DataTypes, ElementType>::addDForce(VecDeriv& df, const VecDeriv& dx,
                                                        Real kFactor) const
{
    const auto& elements = FiniteElement::getElementSequence(*m_topology);

    Deriv nodedForce(sofa::type::NOINIT);

    auto elementStiffnessIt = this->stiffnessMatrices().begin();
    auto rotationMatrixIt = m_rotations.begin();
    for (const auto& element : elements)
    {
        const auto& elementRotation = *rotationMatrixIt++;

        sofa::type::Vec<NumberOfDofsInElement, Real> element_dx;
        for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
        {
            const auto rotated_dx = elementRotation.transposed() * dx[element[i]];
            for (sofa::Size j = 0; j < spatial_dimensions; ++j)
            {
                element_dx[i * spatial_dimensions + j] = rotated_dx[j];
            }
        }

        const auto& stiffnessMatrix = *elementStiffnessIt++;
        const auto dForce = kFactor * stiffnessMatrix * element_dx;

        for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
        {
            dForce.getsub(i * spatial_dimensions, nodedForce);
            df[element[i]] += -elementRotation * nodedForce;
        }
    }
}

template <class DataTypes, class ElementType>
void CorotationalFEM<DataTypes, ElementType>::buildStiffnessMatrix(
    sofa::core::behavior::StiffnessMatrix::Derivative& dfdx) const
{
}

template <class DataTypes, class ElementType>
auto CorotationalFEM<DataTypes, ElementType>::translation(
    const std::array<Coord, NumberOfNodesInElement>& nodes) const -> Coord
{
    // return nodes[0];
    return computeCentroid(nodes);
}

template <class DataTypes, class ElementType>
void CorotationalFEM<DataTypes, ElementType>::applyRotation(ElementStiffness& K, const sofa::type::Quat<Real>& rotation)
{
    sofa::type::Mat<3,3,Real> R(sofa::type::NOINIT);
    rotation.toMatrix(R);

    sofa::type::Mat<3,3,Real> K_ij(sofa::type::NOINIT);
    for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
    {
        for (sofa::Size j = 0; j < NumberOfNodesInElement; ++j)
        {
            K.getsub(i * spatial_dimensions, j * spatial_dimensions, K_ij);
            K_ij = R * K_ij * R.transposed();
            K.setsub(i * spatial_dimensions, j * spatial_dimensions, K_ij);
        }
    }
}

template <class DataTypes, class ElementType>
auto CorotationalFEM<DataTypes, ElementType>::referenceElementCentroid() -> Coord
{
    static const Coord centroid = computeCentroid(FiniteElement::referenceElementNodes);
    return centroid;
}

template <class DataTypes, class ElementType>
auto CorotationalFEM<DataTypes, ElementType>::computeCentroid(
    const std::array<Coord, NumberOfNodesInElement>& nodes) -> Coord
{
    Coord centroid;
    for (const auto node : nodes)
    {
        centroid += node;
    }
    centroid /= static_cast<Real>(FiniteElement::NumberOfNodesInElement);
    return centroid;
}

template <class DataTypes, class ElementType>
auto CorotationalFEM<DataTypes, ElementType>::computeElementRotation(
    const std::array<Coord, NumberOfNodesInElement>& nodesPosition,
    const std::array<Coord, NumberOfNodesInElement>& nodesRestPosition) -> RotationMatrix
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

    const auto H = P.transposed() * Q;

    RotationMatrix rotation;
    sofa::helper::Decompose<Real>::polarDecomposition(H, rotation);
    return rotation;
}

}  // namespace elasticity

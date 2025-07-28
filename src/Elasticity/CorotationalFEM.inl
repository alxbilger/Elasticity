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
        computeElementRotation(elementNodesCoordinates, restElementNodesCoordinates, elementRotation);

        const auto t = translation(elementNodesCoordinates);

        ElementDisplacement displacement(sofa::type::NOINIT);
        for (sofa::Size j = 0; j < NumberOfNodesInElement; ++j)
        {
            const Coord rotatedDisplacement = elementRotation.transposed() * (elementNodesCoordinates[j] - t) - restElementNodesCoordinates[j];
            for (sofa::Size k = 0; k < spatial_dimensions; ++k)
            {
                displacement[j * spatial_dimensions + k] = rotatedDisplacement[k];
            }
        }

        const ElementStiffness& stiffnessMatrix = *elementStiffnessIt++;
        const sofa::type::Vec<NumberOfDofsInElement, Real> elementForce = stiffnessMatrix * displacement;

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
        const auto elementRotation_T = elementRotation.transposed();

        sofa::type::Vec<NumberOfDofsInElement, Real> element_dx(sofa::type::NOINIT);
        for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
        {
            const auto rotated_dx = elementRotation_T * dx[element[i]];
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
    sofa::type::Mat<spatial_dimensions, spatial_dimensions, Real> localMatrix(sofa::type::NOINIT);

    const auto& elements = FiniteElement::getElementSequence(*m_topology);
    auto elementStiffnessIt = this->stiffnessMatrices().begin();
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
                stiffnessMatrix.getsub(spatial_dimensions * n1, spatial_dimensions * n2,
                                       localMatrix);  // extract the submatrix corresponding to the
                                                      // coupling of nodes n1 and n2
                dfdx(element[n1] * spatial_dimensions, element[n2] * spatial_dimensions) +=
                    -elementRotation * localMatrix * elementRotation_T;
            }
        }
    }
}

template <class DataTypes, class ElementType>
void CorotationalFEM<DataTypes, ElementType>::computeVonMisesStress(
    VonMisesStressContainer<Real>& vonMisesStressContainer,
    const VecCoord& position, const VecCoord& restPosition) const
{
    if constexpr (spatial_dimensions > 1)
    {
        const auto& elements = FiniteElement::getElementSequence(*m_topology);

        for (sofa::Size i = 0; i < elements.size(); ++i)
        {
            const auto& element = elements[i];
            const auto& elementRotation = m_rotations[i];

            const std::array<Coord, NumberOfNodesInElement> elementNodesCoordinates = extractNodesVectorFromGlobalVector(element, position);
            const std::array<Coord, NumberOfNodesInElement> restElementNodesCoordinates = extractNodesVectorFromGlobalVector(element, restPosition);

            const auto t = translation(elementNodesCoordinates);

            ElementDisplacement displacement(sofa::type::NOINIT);
            for (sofa::Size j = 0; j < NumberOfNodesInElement; ++j)
            {
                const Coord rotatedDisplacement = elementRotation.transposed() * (elementNodesCoordinates[j] - t) - restElementNodesCoordinates[j];
                for (sofa::Size k = 0; k < spatial_dimensions; ++k)
                {
                    displacement[j * spatial_dimensions + k] = rotatedDisplacement[k];
                }
            }

            for (sofa::Size j = 0; j < NumberOfNodesInElement; ++j)
            {
                const auto& B = this->m_strainDisplacement[i][j];
                const auto strain = B * displacement;
                const auto cauchyStress = this->m_elasticityTensor * strain;
                const auto traceCauchyStress = traceFromVoigtTensor(cauchyStress);

                auto deviatoricStress = cauchyStress;
                for (sofa::Size k = 0; k < spatial_dimensions; ++k)
                    deviatoricStress[k] -= (1./spatial_dimensions) * traceCauchyStress;

                Real vonMisesStressValue = 0;
                for (sofa::Size k = 0; k < spatial_dimensions; ++k)
                    vonMisesStressValue += deviatoricStress[k] * deviatoricStress[k];
                for (sofa::Size k = spatial_dimensions; k < NumberOfIndependentElements; ++k)
                    vonMisesStressValue += 2 * deviatoricStress[k] * deviatoricStress[k];
                vonMisesStressValue = sqrt(static_cast<Real>(spatial_dimensions) / (2. * (spatial_dimensions-1.)) * vonMisesStressValue);
                vonMisesStressContainer.addVonMisesStress(element[j], vonMisesStressValue);
            }
        }
    }
}

template <class DataTypes, class ElementType>
auto CorotationalFEM<DataTypes, ElementType>::translation(
    const std::array<Coord, NumberOfNodesInElement>& nodes) const -> Coord
{
    // return nodes[0];
    return computeCentroid(nodes);
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
void CorotationalFEM<DataTypes, ElementType>::computeElementRotation(
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

    const auto H = P.transposed() * Q;

    sofa::helper::Decompose<Real>::polarDecomposition_stable(H, rotationMatrix);
}

}  // namespace elasticity

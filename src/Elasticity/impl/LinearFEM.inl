#pragma once
#include <Elasticity/impl/ElasticityTensor.h>
#include <Elasticity/impl/ElementStiffnessMatrix.h>
#include <Elasticity/impl/LinearFEM.h>
#include <Elasticity/impl/MatrixTools.h>
#include <Elasticity/impl/VectorTools.h>
#include <sofa/core/behavior/BaseLocalForceFieldMatrix.h>
#include <Elasticity/impl/VecView.h>

namespace elasticity
{

template <class DataTypes, class ElementType, ComputationStrategy strategy>
LinearFEM<DataTypes, ElementType, strategy>::LinearFEM(sofa::core::topology::BaseMeshTopology* topology)
    : m_topology(topology)
{
}

template <class DataTypes, class ElementType, ComputationStrategy strategy>
void LinearFEM<DataTypes, ElementType, strategy>::addForce(VecDeriv& force, const VecCoord& position,
                                                 const VecCoord& restPosition)
{
    if (m_topology == nullptr) return;

    const auto& elements = FiniteElement::getElementSequence(*m_topology);

    for (sofa::Size i = 0; i < elements.size(); ++i)
    {
        const auto& element = elements[i];
        const ElementStiffness& stiffnessMatrix = m_elementStiffness[i];

        const std::array<Coord, NumberOfNodesInElement> elementNodesCoordinates = extractNodesVectorFromGlobalVector(element, position);
        const std::array<Coord, NumberOfNodesInElement> restElementNodesCoordinates = extractNodesVectorFromGlobalVector(element, restPosition);

        const ElementDisplacement displacement = computeElementDisplacement(elementNodesCoordinates, restElementNodesCoordinates);

        sofa::type::Vec<NumberOfDofsInElement, Real> elementForce = stiffnessMatrix * displacement;

        for (sofa::Size j = 0; j < NumberOfNodesInElement; ++j)
        {
            VecView<spatial_dimensions, Real> nodeForce(elementForce, j * spatial_dimensions);
            force[element[j]] += -nodeForce;
        }
    }
}

template <class DataTypes, class ElementType, ComputationStrategy strategy>
void LinearFEM<DataTypes, ElementType, strategy>::addDForce(VecDeriv& df, const VecDeriv& dx, Real kFactor) const
{
    const auto& elements = FiniteElement::getElementSequence(*m_topology);

    auto elementStiffnessIt = stiffnessMatrices().begin();
    for (const auto& element : elements)
    {
        sofa::type::Vec<NumberOfDofsInElement, Real> element_dx(sofa::type::NOINIT);
        for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
        {
            for (sofa::Size j = 0; j < spatial_dimensions; ++j)
            {
                element_dx[i * spatial_dimensions + j] = dx[element[i]][j];
            }
        }

        const auto& stiffnessMatrix = *elementStiffnessIt++;
        auto dForce = (-kFactor) * (stiffnessMatrix * element_dx);

        for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
        {
            VecView<spatial_dimensions, Real> nodedForce(dForce, i * spatial_dimensions);
            df[element[i]] += nodedForce.toVec();
        }
    }
}

template <class DataTypes, class ElementType, ComputationStrategy strategy>
void LinearFEM<DataTypes, ElementType, strategy>::buildStiffnessMatrix(
    sofa::core::behavior::StiffnessMatrix::Derivative& dfdx) const
{
    sofa::type::Mat<spatial_dimensions, spatial_dimensions, Real> localMatrix(sofa::type::NOINIT);

    const auto& elements = FiniteElement::getElementSequence(*m_topology);
    auto elementStiffnessIt = stiffnessMatrices().begin();
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

template <class DataTypes, class ElementType, ComputationStrategy strategy>
void LinearFEM<DataTypes, ElementType, strategy>::precomputeElementStiffness(const VecCoord& restPosition,
                                                                   Real youngModulus,
                                                                   Real poissonRatio)
{
    m_elasticityTensor = makeIsotropicElasticityTensor<DataTypes, strategy>(youngModulus, poissonRatio);

    m_elementStiffness.clear();

    const auto& elements = FiniteElement::getElementSequence(*m_topology);
    m_elementStiffness.reserve(elements.size());

    for (const auto& element : elements)
    {
        const std::array<Coord, NumberOfNodesInElement> nodesCoordinates = extractNodesVectorFromGlobalVector(element, restPosition);
        ElementStiffness K = integrate<DataTypes, ElementType, strategy>(nodesCoordinates, m_elasticityTensor);
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
        const auto nodesCoordinates = extractNodesVectorFromGlobalVector(element, restPosition);
        for (sofa::Size j = 0; j < NumberOfNodesInElement; ++j)
        {
            const ReferenceCoord& x = referenceElementNodes[j];
            // const auto [B, detJ] = computeStrainDisplacement(nodesCoordinates, x);
            // m_strainDisplacement[i][j] = B;
        }
    }
}

template <class DataTypes, class ElementType, ComputationStrategy strategy>
void LinearFEM<DataTypes, ElementType, strategy>::computeVonMisesStress(
    VonMisesStressContainer<Real>& vonMisesStressContainer,
    const VecCoord& position,
    const VecCoord& restPosition) const
{
    if constexpr (spatial_dimensions > 1)
    {
        const auto& elements = FiniteElement::getElementSequence(*m_topology);

        for (sofa::Size i = 0; i < elements.size(); ++i)
        {
            const auto& element = elements[i];

            const std::array<Coord, NumberOfNodesInElement> elementNodesCoordinates = extractNodesVectorFromGlobalVector(element, position);
            const std::array<Coord, NumberOfNodesInElement> restElementNodesCoordinates = extractNodesVectorFromGlobalVector(element, restPosition);

            const ElementDisplacement displacement = computeElementDisplacement(elementNodesCoordinates, restElementNodesCoordinates);

            for (sofa::Size j = 0; j < NumberOfNodesInElement; ++j)
            {
                const auto& B = m_strainDisplacement[i][j];
                const auto strain = B * displacement;
                const auto cauchyStress = m_elasticityTensor * strain;
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

template <class DataTypes, class ElementType, ComputationStrategy strategy>
auto LinearFEM<DataTypes, ElementType, strategy>::computeElementDisplacement(
    const std::array<Coord, NumberOfNodesInElement>& elementNodesCoordinates,
    const std::array<Coord, NumberOfNodesInElement>& restElementNodesCoordinates)
    -> ElementDisplacement
{
    ElementDisplacement displacement;
    for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
    {
        for (sofa::Size j = 0; j < spatial_dimensions; ++j)
        {
            displacement[i * spatial_dimensions + j] = elementNodesCoordinates[i][j] - restElementNodesCoordinates[i][j];
        }
    }
    return displacement;
}

template <class DataTypes, class ElementType, ComputationStrategy strategy>
auto LinearFEM<DataTypes, ElementType, strategy>::stiffnessMatrices() const -> const sofa::type::vector<ElementStiffness>&
{
    return m_elementStiffness;
}

}  // namespace elasticity

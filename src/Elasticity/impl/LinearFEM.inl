#pragma once
#include <Elasticity/impl/LinearFEM.h>
#include <Elasticity/impl/MatrixTools.h>
#include <Elasticity/impl/VectorTools.h>
#include <sofa/core/behavior/BaseLocalForceFieldMatrix.h>

namespace elasticity
{

template <class DataTypes, class ElementType>
LinearFEM<DataTypes, ElementType>::LinearFEM(sofa::core::topology::BaseMeshTopology* topology)
    : m_topology(topology)
{
}

template <class DataTypes, class ElementType>
void LinearFEM<DataTypes, ElementType>::addForce(VecDeriv& force, const VecCoord& position,
                                                 const VecCoord& restPosition)
{
    if (m_topology == nullptr) return;

    const auto& elements = FiniteElement::getElementSequence(*m_topology);

    Deriv nodeForce(sofa::type::NOINIT);

    for (sofa::Size i = 0; i < elements.size(); ++i)
    {
        const auto& element = elements[i];
        const ElementStiffness& stiffnessMatrix = m_elementStiffness[i];

        const std::array<Coord, NumberOfNodesInElement> elementNodesCoordinates = extractNodesVectorFromGlobalVector(element, position);
        const std::array<Coord, NumberOfNodesInElement> restElementNodesCoordinates = extractNodesVectorFromGlobalVector(element, restPosition);

        const ElementDisplacement displacement = computeElementDisplacement(elementNodesCoordinates, restElementNodesCoordinates);

        const sofa::type::Vec<NumberOfDofsInElement, Real> elementForce = stiffnessMatrix * displacement;

        for (sofa::Size j = 0; j < NumberOfNodesInElement; ++j)
        {
            elementForce.getsub(j * spatial_dimensions, nodeForce);
            force[element[j]] += -nodeForce;
        }
    }
}

template <class DataTypes, class ElementType>
void LinearFEM<DataTypes, ElementType>::addDForce(VecDeriv& df, const VecDeriv& dx, Real kFactor) const
{
    const auto& elements = FiniteElement::getElementSequence(*m_topology);

    Deriv nodedForce(sofa::type::NOINIT);

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
        const auto dForce = kFactor * stiffnessMatrix * element_dx;

        for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
        {
            dForce.getsub(i * spatial_dimensions, nodedForce);
            df[element[i]] += -nodedForce;
        }
    }
}

template <class DataTypes, class ElementType>
void LinearFEM<DataTypes, ElementType>::buildStiffnessMatrix(
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

template <class DataTypes, class ElementType>
void LinearFEM<DataTypes, ElementType>::precomputeElementStiffness(const VecCoord& restPosition,
                                                                   Real youngModulus,
                                                                   Real poissonRatio)
{
    m_elasticityTensor = computeElasticityTensor(youngModulus, poissonRatio);

    m_elementStiffness.clear();

    const auto& elements = FiniteElement::getElementSequence(*m_topology);
    m_elementStiffness.reserve(elements.size());

    for (const auto& element : elements)
    {
        ElementStiffness K;
        const auto nodesCoordinates = extractNodesVectorFromGlobalVector(element, restPosition);
        for (const auto& [quadraturePoint, weight] : FiniteElement::quadraturePoints())
        {
            const auto [B, detJ] = computeStrainDisplacement(nodesCoordinates, quadraturePoint);

            K += (weight * detJ) * B.transposed() * m_elasticityTensor * B;
        }
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
            const auto [B, detJ] = computeStrainDisplacement(nodesCoordinates, x);
            m_strainDisplacement[i][j] = B;
        }
    }
}

template <class DataTypes, class ElementType>
void LinearFEM<DataTypes, ElementType>::computeVonMisesStress(
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

template <class DataTypes, class ElementType>
typename LinearFEM<DataTypes, ElementType>::ElasticityTensor
LinearFEM<DataTypes, ElementType>::computeElasticityTensor(Real youngModulus, Real poissonRatio)
{
    static constexpr auto volumetricTensor = []()
    {
        sofa::type::Mat<NumberOfIndependentElements, NumberOfIndependentElements, Real> I_vol;
        for (sofa::Size i = 0; i < spatial_dimensions; ++i)
        {
            for (sofa::Size j = 0; j < spatial_dimensions; ++j)
            {
                I_vol[i][j] = 1;
            }
        }
        return I_vol;
    }();

    static constexpr auto deviatoricTensor = []()
    {
        sofa::type::Mat<NumberOfIndependentElements, NumberOfIndependentElements, Real> I_dev;
        for (sofa::Size i = 0; i < spatial_dimensions; ++i)
        {
            I_dev[i][i] = 1;
        }
        for (sofa::Size i = spatial_dimensions; i < NumberOfIndependentElements; ++i)
        {
            I_dev[i][i] = 0.5;
        }
        return I_dev;
    }();

    const auto mu = youngModulus / (2 * (1 + poissonRatio));
    const auto lambda = youngModulus * poissonRatio /
                        ((1 + poissonRatio) * (1 - (spatial_dimensions - 1) * poissonRatio));

    return lambda * volumetricTensor + 2 * mu * deviatoricTensor;
}

template <class DataTypes, class ElementType>
auto LinearFEM<DataTypes, ElementType>::buildStrainDisplacement(
    const sofa::type::Mat<NumberOfNodesInElement, spatial_dimensions, Real> gradientShapeFunctions)
-> StrainDisplacement
{
    StrainDisplacement B;
    for (sofa::Size ne = 0; ne < NumberOfNodesInElement; ++ne)
    {
        for (sofa::Size i = 0; i < spatial_dimensions; ++i)
        {
            B[i][ne * spatial_dimensions + i] = gradientShapeFunctions[ne][i];
        }

        auto row = spatial_dimensions;
        for (sofa::Size i = 0; i < spatial_dimensions; ++i)
        {
            for (sofa::Size j = i + 1; j < spatial_dimensions; ++j)
            {
                B[row][ne * spatial_dimensions + i] = gradientShapeFunctions[ne][j];
                B[row][ne * spatial_dimensions + j] = gradientShapeFunctions[ne][i];
                ++row;
            }
        }
    }

    return B;
}

template <class DataTypes, class ElementType>
auto LinearFEM<DataTypes, ElementType>::computeElementDisplacement(
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

template <class DataTypes, class ElementType>
auto LinearFEM<DataTypes, ElementType>::stiffnessMatrices() const -> const sofa::type::vector<ElementStiffness>&
{
    return m_elementStiffness;
}

template <class DataTypes, class ElementType>
auto LinearFEM<DataTypes, ElementType>::computeStrainDisplacement(
    const std::array<Coord, NumberOfNodesInElement>& nodesCoordinates,
    const ReferenceCoord& evaluationPoint) -> std::pair<StrainDisplacement, Real>
{
    // gradient of shape functions in the reference element evaluated at the quadrature point
    const sofa::type::Mat<NumberOfNodesInElement, ElementDimension, Real> dN_dq_ref =
        FiniteElement::gradientShapeFunctions(evaluationPoint);

    // jacobian of the mapping from the reference space to the physical space, evaluated at the
    // quadrature point
    sofa::type::Mat<spatial_dimensions, ElementDimension, Real> jacobian;
    for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
        jacobian += sofa::type::dyad(nodesCoordinates[i], dN_dq_ref[i]);

    const auto detJ = elasticity::determinant(jacobian);
    const sofa::type::Mat<ElementDimension, spatial_dimensions, Real> J_inv =
        elasticity::inverse(jacobian);

    // gradient of the shape functions in the physical element evaluated at the quadrature point
    sofa::type::Mat<NumberOfNodesInElement, spatial_dimensions, Real> dN_dq(sofa::type::NOINIT);
    for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
        dN_dq[i] = J_inv.transposed() * dN_dq_ref[i];

    return {buildStrainDisplacement(dN_dq), detJ};
}

}  // namespace elasticity

#pragma once
#include <Elasticity/LinearFEM.h>
#include <Elasticity/MatrixTools.h>
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

    updateStiffnessMatrices(position, restPosition);
    auto elementStiffnessIt = stiffnessMatrices().begin();

    for (const auto& element : elements)
    {
        std::array<Coord, NumberOfNodesInElement> elementNodesCoordinates, restElementNodesCoordinates;
        for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
        {
            elementNodesCoordinates[i] = position[element[i]];
            restElementNodesCoordinates[i] = restPosition[element[i]];
        }

        const ElementDisplacement displacement = computeElementDisplacement(elementNodesCoordinates, restElementNodesCoordinates);

        const ElementStiffness& stiffnessMatrix = *elementStiffnessIt++;
        const sofa::type::Vec<NumberOfDofsInElement, Real> elementForce = stiffnessMatrix * displacement;

        for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
        {
            elementForce.getsub(i * spatial_dimensions, nodeForce);
            force[element[i]] += -nodeForce;
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
        sofa::type::Vec<NumberOfDofsInElement, Real> element_dx;
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
void LinearFEM<DataTypes, ElementType>::precomputeElementStiffness(const VecCoord& restPosition, Real youngModulus, Real poissonRatio)
{
    const ElasticityTensor C = computeElasticityTensor(youngModulus, poissonRatio);

    m_elementStiffness.clear();

    const auto& elements = FiniteElement::getElementSequence(*m_topology);
    m_elementStiffness.reserve(elements.size());

    for (const auto& element : elements)
    {
        // matrix where the i-th column is the i-th node coordinates in the element
        const sofa::type::Mat<spatial_dimensions, NumberOfNodesInElement, Real> X_element =
            nodesMatrix(element, restPosition);

        ElementStiffness K;
        for (const auto& [quadraturePoint, weight] : FiniteElement::quadraturePoints())
        {
            // gradient of shape functions in the reference element evaluated at the quadrature
            // point
            const sofa::type::Mat<NumberOfNodesInElement, ElementDimension, Real> dN_dq_ref =
                FiniteElement::gradientShapeFunctions(quadraturePoint);

            // jacobian of the mapping from the reference space to the physical space, evaluated at
            // the quadrature point
            const sofa::type::Mat<spatial_dimensions, ElementDimension, Real> jacobian =
                X_element * dN_dq_ref;

            const auto detJ = elasticity::determinant(jacobian);
            const sofa::type::Mat<ElementDimension, spatial_dimensions, Real> J_inv =
                elasticity::inverse(jacobian);

            // gradient of the shape functions in the physical element evaluated at the quadrature
            // point
            const sofa::type::Mat<NumberOfNodesInElement, spatial_dimensions, Real> dN_dq =
                dN_dq_ref * J_inv;

            const auto B = buildStrainDisplacement(dN_dq);

            K += (weight * detJ) * B.transposed() * C * B;
        }
        m_elementStiffness.push_back(K);
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
            displacement[i * spatial_dimensions + j] =
                elementNodesCoordinates[i][j] - restElementNodesCoordinates[i][j];
        }
    }
    return displacement;
}

template <class DataTypes, class ElementType>
void LinearFEM<DataTypes, ElementType>::updateStiffnessMatrices(const VecCoord& positions,
                                                                const VecCoord& restPositions)
{
    SOFA_UNUSED(positions);
    SOFA_UNUSED(restPositions);
}

template <class DataTypes, class ElementType>
auto LinearFEM<DataTypes, ElementType>::stiffnessMatrices() const -> const sofa::type::vector<ElementStiffness>&
{
    return m_elementStiffness;
}

}  // namespace elasticity

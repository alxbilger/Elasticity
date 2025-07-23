#pragma once
#include <Elasticity/CorotationalFEM.h>
#include <Elasticity/MatrixTools.h>
#include <Elasticity/VectorTools.h>
#include <sofa/helper/decompose.h>

namespace elasticity
{

template <class DataTypes, class ElementType>
void CorotationalFEM<DataTypes, ElementType>::updateStiffnessMatrices(const VecCoord& positions, const VecCoord& restPositions)
{
    SOFA_UNUSED(restPositions);

    const auto& elements = FiniteElement::getElementSequence(*this->m_topology);
    assert(this->m_elementStiffness.size() == elements.size());

    m_rotatedStiffness.resize(elements.size());

    for (sofa::Size i = 0; i < elements.size(); ++i)
    {
        const auto& element = elements[i];

        ElementStiffness& K = m_rotatedStiffness[i];
        K = this->m_elementStiffness[i];

        std::array<Coord, NumberOfNodesInElement> nodes = extractNodesVectorFromGlobalVector(element, positions);
        const auto centroid = computeCentroid(nodes);

        std::array<Coord, NumberOfNodesInElement> restnodes = extractNodesVectorFromGlobalVector(element, restPositions);
        const auto restcentroid = computeCentroid(restnodes);

        sofa::type::Mat<NumberOfNodesInElement, spatial_dimensions, Real> P(sofa::type::NOINIT);
        sofa::type::Mat<NumberOfNodesInElement, spatial_dimensions, Real> Q(sofa::type::NOINIT);

        for (sofa::Size j = 0; j < NumberOfNodesInElement; ++j)
        {
            P[j] = nodes[j] - centroid;
            Q[j] = restnodes[j] - restcentroid;
        }

        const auto H = P.transposed() * Q;



        //
        // const std::array<Coord, NumberOfNodesInElement> nodes = extractNodesVectorFromGlobalVector(element, positions);
        // const auto centroid = computeCentroid(nodes);
        //
        // sofa::type::Mat<NumberOfNodesInElement, spatial_dimensions, Real> X_T(sofa::type::NOINIT);
        // for (sofa::Size j = 0; j < NumberOfNodesInElement; ++j)
        // {
        //     X_T[j] = nodes[j] - centroid;
        // }
        //
        // // matrix where the i-th column is the i-th node coordinates in the element
        // const sofa::type::Mat<spatial_dimensions, NumberOfNodesInElement, Real> X_element = X_T.transposed();
        //
        // const DeformationGradient F = deformationGradient(X_element, m_restJacobians[i]);
        //
        // // ::elasticity::extractRotation(F, m_rotations[i], 1000);

        sofa::type::Mat<3,3,Real> R_0_1;
        sofa::helper::Decompose<Real>::polarDecomposition( H, R_0_1 );
        m_rotations[i].fromMatrix(R_0_1);

        applyRotation(K, m_rotations[i]);
    }
}

template <class DataTypes, class ElementType>
auto CorotationalFEM<DataTypes, ElementType>::stiffnessMatrices() const
    -> const sofa::type::vector<ElementStiffness>&
{
    return m_rotatedStiffness;
}

template <class DataTypes, class ElementType>
auto CorotationalFEM<DataTypes, ElementType>::deformationGradient(
    const sofa::type::Mat<spatial_dimensions, NumberOfNodesInElement, Real>& nodesMatrix,
    const sofa::type::Mat<spatial_dimensions, ElementDimension, Real>& inverseJacobian) -> DeformationGradient
{
    // gradient of shape functions in the reference element evaluated at the centroid of the element
    static const sofa::type::Mat<NumberOfNodesInElement, ElementDimension, Real> dN_dq_ref =
        FiniteElement::gradientShapeFunctions(referenceElementCentroid());

    // jacobian of the mapping from the reference space to the physical space, evaluated at the
    // centroid of the physical element
    const sofa::type::Mat<spatial_dimensions, ElementDimension, Real> jacobian = nodesMatrix * dN_dq_ref;

    return jacobian * inverseJacobian;
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
void CorotationalFEM<DataTypes, ElementType>::computeRestJacobians(const VecCoord& restPosition)
{
    m_restJacobians.clear();

    const auto& elements = FiniteElement::getElementSequence(*this->m_topology);
    m_restJacobians.reserve(elements.size());

    for (const auto& element : elements)
    {
        const std::array<Coord, NumberOfNodesInElement> nodes = extractNodesVectorFromGlobalVector(element, restPosition);
        const auto centroid = computeCentroid(nodes);

        sofa::type::Mat<NumberOfNodesInElement, spatial_dimensions, Real> X_T(sofa::type::NOINIT);
        for (sofa::Size j = 0; j < NumberOfNodesInElement; ++j)
        {
            X_T[j] = nodes[j] - centroid;
        }

        // matrix where the i-th column is the i-th node coordinates in the element
        const sofa::type::Mat<spatial_dimensions, NumberOfNodesInElement, Real> X_element = X_T.transposed();

        // gradient of shape functions in the reference element evaluated at the centroid of the
        // element
        static const sofa::type::Mat<NumberOfNodesInElement, ElementDimension, Real> dN_dq_ref =
            FiniteElement::gradientShapeFunctions(referenceElementCentroid());

        // jacobian of the mapping from the reference space to the physical space, evaluated at the
        // centroid of the rest element
        const sofa::type::Mat<spatial_dimensions, ElementDimension, Real> jacobian =
            X_element * dN_dq_ref;

        m_restJacobians.push_back(jacobian.inverted());
    }
}

template <class DataTypes, class ElementType>
void CorotationalFEM<DataTypes, ElementType>::precomputeElementStiffness(
    const VecCoord& restPosition, Real youngModulus, Real poissonRatio)
{
    LinearFEM<DataTypes, ElementType>::precomputeElementStiffness(restPosition, youngModulus, poissonRatio);

    computeRestJacobians(restPosition);

    const auto& elements = FiniteElement::getElementSequence(*this->m_topology);
    m_rotations.resize(elements.size());
}

}  // namespace elasticity

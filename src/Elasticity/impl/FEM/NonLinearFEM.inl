#pragma once
#include <Elasticity/impl/ElasticityTensor.h>
#include <Elasticity/impl/FEM/NonLinearFEM.h>
#include <Elasticity/impl/MatrixTools.h>
#include <Elasticity/impl/VectorTools.h>

namespace elasticity
{

template <class DataTypes, class ElementType>
NonLinearFEM<DataTypes, ElementType>::NonLinearFEM(sofa::core::topology::BaseMeshTopology* topology)
    : m_topology(topology)
{
}

template <class DataTypes, class ElementType>
void NonLinearFEM<DataTypes, ElementType>::addForce(VecDeriv& force, const VecCoord& position,
                                                    const VecCoord& restPosition)
{
    if (m_topology == nullptr) return;

    const auto& elements = FiniteElement::getElementSequence(*m_topology);

    for (const auto& element : elements)
    {
        const std::array<Coord, NumberOfNodesInElement> elementNodesCoordinates = extractNodesVectorFromGlobalVector(element, position);
        const std::array<Coord, NumberOfNodesInElement> elementNodesRestCoordinates = extractNodesVectorFromGlobalVector(element, restPosition);

        for (const auto& [quadraturePoint, weight] : FiniteElement::quadraturePoints())
        {
            // gradient of shape functions in the reference element evaluated at the quadrature point
            const sofa::type::Mat<NumberOfNodesInElement, ElementDimension, Real> dN_dq_ref =
                FiniteElement::gradientShapeFunctions(quadraturePoint);

            // jacobian of the mapping from the reference space to the physical space, evaluated at the
            // quadrature point
            sofa::type::Mat<spatial_dimensions, ElementDimension, Real> J_q;
            for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
                J_q += sofa::type::dyad(elementNodesCoordinates[i], dN_dq_ref[i]);

            // jacobian of the mapping from the reference space to the rest physical space, evaluated at the
            // quadrature point
            sofa::type::Mat<spatial_dimensions, ElementDimension, Real> J_Q;
            for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
                J_Q += sofa::type::dyad(elementNodesRestCoordinates[i], dN_dq_ref[i]);

            const auto detJ_q = elasticity::determinant(J_q);
            const auto detJ_Q = elasticity::determinant(J_Q);

            const sofa::type::Mat<ElementDimension, spatial_dimensions, Real> J_q_inv =
                elasticity::inverse(J_q);
            const sofa::type::Mat<ElementDimension, spatial_dimensions, Real> J_Q_inv =
                elasticity::inverse(J_Q);

            // gradient of the shape functions in the physical element evaluated at the quadrature point
            sofa::type::Mat<NumberOfNodesInElement, spatial_dimensions, Real> dN_dq(sofa::type::NOINIT);
            for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
                dN_dq[i] = J_q_inv.transposed() * dN_dq_ref[i];

            // gradient of the shape functions in the physical element evaluated at the quadrature point
            sofa::type::Mat<NumberOfNodesInElement, spatial_dimensions, Real> dN_dQ(sofa::type::NOINIT);
            for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
                dN_dQ[i] = J_Q_inv.transposed() * dN_dq_ref[i];

            DeformationGradient F = J_q * elasticity::inverse(J_Q);

            const auto detF = elasticity::determinant(F);
            if (detF < 0)
            {
                msg_error("FEM") << "Element inversion detected (detF = " << detF << " < 0, " <<
                    " detJ_q = " << detJ_q << ", detJ_Q = " << elasticity::determinant(J_Q) << ")";
            }

            // Right Cauchy-Green deformation tensor
            const auto C = F.transposed() * F;

            static const auto& I = sofa::type::Mat<spatial_dimensions, spatial_dimensions, Real>::Identity();

            // Green-Lagrangian strain tensor
            const auto E = 0.5 * (C - I);

            const auto [lambda, mu] = elasticity::toLameParameters<DataTypes>(10000., 0.45);

            // First Piola-Kirchhoff stress tensor
            const auto S = lambda * sofa::type::trace(E) * I + 2 * mu * E;

            // Second Piola-Kirchhoff stress tensor
            const auto P = F * S;

            for (sofa::Index i = 0; i < NumberOfNodesInElement; ++i)
            {
                force[element[i]] += (-detJ_Q * weight) * P * dN_dQ[i];
            }
        }
    }
}

template <class DataTypes, class ElementType>
void NonLinearFEM<DataTypes, ElementType>::addDForce(VecDeriv& df, const VecDeriv& dx, Real kFactor) const
{
    if (m_topology == nullptr) return;

    const auto& elements = FiniteElement::getElementSequence(*m_topology);

    for (const auto& element : elements)
    {
    }
}

template <class DataTypes, class ElementType>
void NonLinearFEM<DataTypes, ElementType>::buildStiffnessMatrix(
    sofa::core::behavior::StiffnessMatrix::Derivative& dfdx) const
{
    if (m_topology == nullptr) return;

    const auto& elements = FiniteElement::getElementSequence(*m_topology);

    for (const auto& element : elements)
    {
    }
}

}  // namespace elasticity

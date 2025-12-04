#pragma once

#include <Elasticity/component/BaseElementLinearFEMForceField.h>
#include <Elasticity/impl/LameParameters.h>
#include <Elasticity/impl/VectorTools.h>

namespace elasticity
{

template <class DataTypes, class ElementType>
void BaseElementLinearFEMForceField<DataTypes, ElementType>::init()
{
    LinearMechanicalParametersComponent<DataTypes>::init();
    sofa::core::behavior::ForceField<DataTypes>::init();
    TopologyAccessor::init();

    if (!this->isComponentStateInvalid())
    {
        this->precomputeElementStiffness();
    }
}

template <class DataTypes, class ElementType>
BaseElementLinearFEMForceField<DataTypes, ElementType>::BaseElementLinearFEMForceField()
{
    this->addUpdateCallback("precomputeStiffness", {&this->d_youngModulus, &this->d_poissonRatio},
    [this](const sofa::core::DataTracker& )
    {
        precomputeElementStiffness();
        return this->getComponentState();
    }, {});
}

template <class DataTypes, class ElementType>
void BaseElementLinearFEMForceField<DataTypes, ElementType>::precomputeElementStiffness()
{
    if (!l_topology)
        return;

    if (this->isComponentStateInvalid())
        return;

    const auto youngModulus = this->d_youngModulus.getValue();
    const auto poissonRatio = this->d_poissonRatio.getValue();
    const auto [mu, lambda] = elasticity::toLameParameters<sofa::defaulttype::Vec3Types>(youngModulus, poissonRatio);

    auto restPositionAccessor = this->mstate->readRestPositions();

    m_elasticityTensor = makeIsotropicElasticityTensor<DataTypes>(mu, lambda);

    m_elementStiffness.clear();

    const auto& elements = trait::FiniteElement::getElementSequence(*l_topology);
    m_elementStiffness.reserve(elements.size());

    for (const auto& element : elements)
    {
        const std::array<sofa::Coord_t<DataTypes>, trait::NumberOfNodesInElement> nodesCoordinates = extractNodesVectorFromGlobalVector(element, restPositionAccessor.ref());
        ElementStiffness K = integrate<DataTypes, ElementType, trait::matrixVectorProductType>(nodesCoordinates, m_elasticityTensor);
        m_elementStiffness.push_back(K);
    }

    // precompute strain-displacement tensors at the nodes positions
    m_strainDisplacement.clear();
    m_strainDisplacement.resize(elements.size());

    static const std::array<typename trait::ReferenceCoord, trait::NumberOfNodesInElement>& referenceElementNodes =
        trait::FiniteElement::referenceElementNodes;

    for (sofa::Size i = 0; i < elements.size(); ++i)
    {
        const auto& element = elements[i];
        const auto nodesCoordinates = extractNodesVectorFromGlobalVector(element, restPositionAccessor.ref());
        for (sofa::Size j = 0; j < trait::NumberOfNodesInElement; ++j)
        {
            const auto& x = referenceElementNodes[j];

            // gradient of shape functions in the reference element evaluated at the quadrature point
            const sofa::type::Mat<trait::NumberOfNodesInElement, trait::TopologicalDimension, sofa::Real_t<DataTypes>> dN_dq_ref =
                trait::FiniteElement::gradientShapeFunctions(x);

            // jacobian of the mapping from the reference space to the physical space, evaluated at the
            // quadrature point
            sofa::type::Mat<trait::spatial_dimensions, trait::TopologicalDimension, sofa::Real_t<DataTypes>> jacobian;
            for (sofa::Size n = 0; n < trait::NumberOfNodesInElement; ++n)
                jacobian += sofa::type::dyad(nodesCoordinates[n], dN_dq_ref[n]);

            const sofa::type::Mat<trait::TopologicalDimension, trait::spatial_dimensions, sofa::Real_t<DataTypes>> J_inv =
                elasticity::inverse(jacobian);

            // gradient of the shape functions in the physical element evaluated at the quadrature point
            sofa::type::Mat<trait::NumberOfNodesInElement, trait::spatial_dimensions, sofa::Real_t<DataTypes>> dN_dq(sofa::type::NOINIT);
            for (sofa::Size n = 0; n < trait::NumberOfNodesInElement; ++n)
                dN_dq[n] = J_inv.transposed() * dN_dq_ref[n];

            const auto B = makeStrainDisplacement<DataTypes, ElementType>(dN_dq);

            m_strainDisplacement[i][j] = B;
        }
    }
}

}

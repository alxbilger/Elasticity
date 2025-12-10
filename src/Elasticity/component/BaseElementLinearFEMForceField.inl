#pragma once

#include <Elasticity/component/BaseElementLinearFEMForceField.h>
#include <Elasticity/impl/LameParameters.h>
#include <Elasticity/impl/VectorTools.h>

#include <execution>
#include <ranges>

#include "FEMForceField.h"

namespace elasticity
{

template <class DataTypes, class ElementType>
void BaseElementLinearFEMForceField<DataTypes, ElementType>::init()
{
    LinearMechanicalParametersComponent<DataTypes>::init();
    TopologyAccessor::init();
    sofa::core::behavior::SingleStateAccessor<DataTypes>::init();

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

    const auto& elements = trait::FiniteElement::getElementSequence(*l_topology);
    m_elementStiffness.resize(elements.size());

    SCOPED_TIMER("precomputeStiffness");
    std::ranges::iota_view indices {static_cast<decltype(elements.size())>(0ul), elements.size()};
    std::for_each(std::execution::par, indices.begin(), indices.end(),
        [&](const auto elementId)
        {
            const auto& element = elements[elementId];
            const std::array<sofa::Coord_t<DataTypes>, trait::NumberOfNodesInElement> nodesCoordinates = extractNodesVectorFromGlobalVector(element, restPositionAccessor.ref());
            m_elementStiffness[elementId] = integrate<DataTypes, ElementType, trait::matrixVectorProductType>(nodesCoordinates, m_elasticityTensor);
        });
}

}

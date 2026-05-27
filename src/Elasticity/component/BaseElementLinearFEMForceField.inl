#pragma once

#include <Elasticity/component/BaseElementLinearFEMForceField.h>
#include <sofa/component/solidmechanics/fem/elastic/BaseLinearElasticityFEMForceField.inl>
#include <sofa/component/solidmechanics/fem/elastic/impl/LameParameters.h>
#include <sofa/component/solidmechanics/fem/elastic/impl/VectorTools.h>

#include <ranges>

#include <Elasticity/component/FEMForceField.h>

namespace elasticity
{

template <class DataTypes, class ElementType>
BaseElementLinearFEMForceField<DataTypes, ElementType>::BaseElementLinearFEMForceField()
    : d_elementStiffness(initData(&d_elementStiffness, "elementStiffness", "List of stiffness matrices per element"))
{
    this->addUpdateCallback("precomputeStiffness", {&this->d_youngModulus, &this->d_poissonRatio},
    [this](const sofa::core::DataTracker& )
    {
        precomputeElementStiffness();
        return this->getComponentState();
    }, {});
}

template <class DataTypes, class ElementType>
void BaseElementLinearFEMForceField<DataTypes, ElementType>::init()
{
    sofa::component::solidmechanics::fem::elastic::BaseLinearElasticityFEMForceField<DataTypes>::init();

    if (!this->isComponentStateInvalid())
    {
        this->precomputeElementStiffness();
    }
}

template <class DataTypes, class ElementType>
void BaseElementLinearFEMForceField<DataTypes, ElementType>::precomputeElementStiffness()
{
    if (!this->l_topology)
        return;

    if (this->isComponentStateInvalid())
        return;

    if (!this->mstate)
        return;

    const auto youngModulusAccessor = sofa::helper::ReadAccessor(this->d_youngModulus);
    const auto poissonRatioAccessor = sofa::helper::ReadAccessor(this->d_poissonRatio);

    auto restPositionAccessor = this->mstate->readRestPositions();

    const auto& elements = trait::FiniteElement::getElementSequence(*this->l_topology);

    auto elementStiffness = sofa::helper::getWriteOnlyAccessor(d_elementStiffness);
    elementStiffness.resize(elements.size());

    SCOPED_TIMER("precomputeStiffness");
    sofa::helper::IotaView indices {static_cast<decltype(elements.size())>(0ul), elements.size()};
    std::for_each(indices.begin(), indices.end(),
        [&](const auto elementId)
        {
            const auto& element = elements[elementId];

            const auto youngModulus = this->getYoungModulusInElement(elementId);
            const auto poissonRatio = this->getPoissonRatioInElement(elementId);

            using Real = sofa::Real_t<DataTypes>;

            sofa::component::solidmechanics::fem::elastic::LameLambda<Real> lambda { 0 };
            sofa::component::solidmechanics::fem::elastic::LameMu<Real> mu { 0 };

            sofa::component::solidmechanics::fem::elastic::toLameParameters<DataTypes::spatial_dimensions, Real>(
                sofa::component::solidmechanics::fem::elastic::YoungModulus<Real>(youngModulus),
                sofa::component::solidmechanics::fem::elastic::PoissonRatio<Real>(poissonRatio),
                lambda, mu);

            const auto elasticityTensor = sofa::component::solidmechanics::fem::elastic::makeIsotropicElasticityTensor<DataTypes::spatial_dimensions>(mu, lambda);

            const std::array<sofa::Coord_t<DataTypes>, trait::NumberOfNodesInElement> nodesCoordinates = sofa::component::solidmechanics::fem::elastic::extractNodesVectorFromGlobalVector(element, restPositionAccessor.ref());
            elementStiffness[elementId] = sofa::component::solidmechanics::fem::elastic::integrate<DataTypes, ElementType, trait::matrixVectorProductType>(nodesCoordinates, elasticityTensor);
        });
}

}

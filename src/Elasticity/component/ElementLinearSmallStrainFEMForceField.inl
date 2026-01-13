#pragma once
#include <Elasticity/component/ElementLinearSmallStrainFEMForceField.h>
#include <Elasticity/component/BaseElementLinearFEMForceField.inl>
#include <Elasticity/component/FEMForceField.inl>
#include <Elasticity/impl/LameParameters.h>
#include <Elasticity/impl/VecView.h>
#include <sofa/core/behavior/ForceField.inl>
#include <ranges>
#include <execution>

namespace elasticity
{

template <class DataTypes, class ElementType>
void ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::init()
{
    BaseElementLinearFEMForceField<DataTypes, ElementType>::init();
    FEMForceField<DataTypes, ElementType>::init();

    if (!this->isComponentStateInvalid())
    {
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Valid);
    }
}


template <class DataTypes, class ElementType>
void ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::computeElementsForces(
    const sofa::simulation::Range<std::size_t>& range,
    const sofa::core::MechanicalParams* mparams,
    sofa::type::vector<ElementForce>& elementForces,
    const sofa::VecCoord_t<DataTypes>& nodePositions)
{
    const auto& elements = trait::FiniteElement::getElementSequence(*this->l_topology);
    auto restPositionAccessor = this->mstate->readRestPositions();

    for (std::size_t elementId = range.start; elementId < range.end; ++elementId)
    {
        const auto& element = elements[elementId];
        const auto& stiffnessMatrix = this->m_elementStiffness[elementId];

        typename trait::ElementDisplacement displacement{ sofa::type::NOINIT };

        for (sofa::Size j = 0; j < trait::NumberOfNodesInElement; ++j)
        {
            const auto nodeId = element[j];
            for (sofa::Size k = 0; k < trait::spatial_dimensions; ++k)
            {
                displacement[j * trait::spatial_dimensions + k] = nodePositions[nodeId][k] - restPositionAccessor[nodeId][k];
            }
        }

        elementForces[elementId] = stiffnessMatrix * displacement;
    }
}

template <class DataTypes, class ElementType>
void ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::computeElementsForcesDeriv(
    const sofa::simulation::Range<std::size_t>& range,
    const sofa::core::MechanicalParams* mparams,
    sofa::type::vector<ElementForce>& elementForcesDeriv,
    const sofa::VecDeriv_t<DataTypes>& nodeDx,
    sofa::Real_t<DataTypes> kFactor)
{
    const auto& elements = trait::FiniteElement::getElementSequence(*this->l_topology);

    for (std::size_t elementId = range.start; elementId < range.end; ++elementId)
    {
        const auto& element = elements[elementId];

        sofa::type::Vec<trait::NumberOfDofsInElement, sofa::Real_t<DataTypes>> element_dx(sofa::type::NOINIT);
        for (sofa::Size i = 0; i < trait::NumberOfNodesInElement; ++i)
        {
            VecView<trait::spatial_dimensions, sofa::Real_t<DataTypes>> node_dx(element_dx, i * trait::spatial_dimensions);
            node_dx = nodeDx[element[i]];
        }

        const auto& stiffnessMatrix = this->m_elementStiffness[elementId];
        elementForcesDeriv[elementId] = kFactor * (stiffnessMatrix * element_dx);
    }
}

template <class DataTypes, class ElementType>
void ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::buildStiffnessMatrix(
    sofa::core::behavior::StiffnessMatrix* matrix)
{
    if (this->isComponentStateInvalid())
        return;

    auto dfdx = matrix->getForceDerivativeIn(this->sofa::core::behavior::ForceField<DataTypes>::mstate)
        .withRespectToPositionsIn(this->sofa::core::behavior::ForceField<DataTypes>::mstate);

    sofa::type::Mat<trait::spatial_dimensions, trait::spatial_dimensions, sofa::Real_t<DataTypes>> localMatrix(sofa::type::NOINIT);

    const auto& elements = trait::FiniteElement::getElementSequence(*this->l_topology);

    if (this->m_elementStiffness.size() < elements.size())
    {
        return;
    }

    auto elementStiffnessIt = this->m_elementStiffness.begin();
    for (const auto& element : elements)
    {
        const auto& stiffnessMatrix = *elementStiffnessIt++;

        for (sofa::Index n1 = 0; n1 < trait::NumberOfNodesInElement; ++n1)
        {
            for (sofa::Index n2 = 0; n2 < trait::NumberOfNodesInElement; ++n2)
            {
                stiffnessMatrix.getAssembledMatrix().getsub(trait::spatial_dimensions * n1, trait::spatial_dimensions * n2, localMatrix); //extract the submatrix corresponding to the coupling of nodes n1 and n2
                dfdx(element[n1] * trait::spatial_dimensions, element[n2] * trait::spatial_dimensions) += -localMatrix;
            }
        }
    }
}

template <class DataTypes, class ElementType>
SReal ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::getPotentialEnergy(
    const sofa::core::MechanicalParams*,
    const sofa::DataVecCoord_t<DataTypes>& x) const
{
    return 0;
}

template <class DataTypes, class ElementType>
void ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::addKToMatrix(
    sofa::linearalgebra::BaseMatrix* matrix, SReal kFact, unsigned& offset)
{
    if (this->isComponentStateInvalid())
        return;

    using LocalMatType = sofa::type::Mat<trait::spatial_dimensions, trait::spatial_dimensions, sofa::Real_t<DataTypes>>;
    LocalMatType localMatrix{sofa::type::NOINIT};

    const auto& elements = trait::FiniteElement::getElementSequence(*this->l_topology);
    auto elementStiffnessIt = this->m_elementStiffness.begin();
    for (const auto& element : elements)
    {
        const auto& stiffnessMatrix = *elementStiffnessIt++;

        for (sofa::Index n1 = 0; n1 < trait::NumberOfNodesInElement; ++n1)
        {
            for (sofa::Index n2 = 0; n2 < trait::NumberOfNodesInElement; ++n2)
            {
                stiffnessMatrix.getAssembledMatrix().getsub(trait::spatial_dimensions * n1, trait::spatial_dimensions * n2, localMatrix); //extract the submatrix corresponding to the coupling of nodes n1 and n2

                const auto value = (-static_cast<sofa::Real_t<DataTypes>>(kFact)) * static_cast<ScalarOrMatrix<LocalMatType>>(localMatrix);
                matrix->add(
                   offset + element[n1] * trait::spatial_dimensions,
                   offset + element[n2] * trait::spatial_dimensions, value);
            }
        }
    }
}

}  // namespace elasticity

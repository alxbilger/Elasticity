#pragma once
#include <Elasticity/component/ElementLinearSmallStrainFEMForceField.h>
#include <Elasticity/component/BaseElementLinearFEMForceField.inl>
#include <Elasticity/impl/LameParameters.h>
#include <Elasticity/impl/VecView.h>
#include <Elasticity/impl/VectorTools.h>
#include <sofa/core/behavior/ForceField.inl>
#include <ranges>
#include <execution>

namespace elasticity
{

constexpr std::string_view sequencedComputeStrategy = "sequenced";
constexpr std::string_view unsequencedComputeStrategy = "unsequenced";
constexpr std::string_view parallelComputeStrategy = "parallel";
constexpr std::string_view parallelUnsequencedComputeStrategy = "parallel_unsequenced";


template <class DataTypes, class ElementType>
ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::ElementLinearSmallStrainFEMForceField()
    : d_computeStrategy(initData(&d_computeStrategy, "computeStrategy", "The compute strategy used to compute the forces"))
{
    sofa::helper::OptionsGroup computeStrategyOptions{
        sequencedComputeStrategy,
        unsequencedComputeStrategy,
        parallelComputeStrategy,
        parallelUnsequencedComputeStrategy
    };
    computeStrategyOptions.setSelectedItem(std::string(parallelComputeStrategy));
    d_computeStrategy.setValue(computeStrategyOptions);

    this->addUpdateCallback("selectStrategy", {&this->d_computeStrategy},
    [this](const sofa::core::DataTracker& )
    {
        selectStrategy();
        return this->getComponentState();
    }, {});
}

template <class DataTypes, class ElementType>
void ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::init()
{
    BaseElementLinearFEMForceField<DataTypes, ElementType>::init();
    sofa::core::behavior::ForceField<DataTypes>::init();

    if (!this->isComponentStateInvalid())
    {
        selectStrategy();
    }

    if (!this->isComponentStateInvalid())
    {
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Valid);
    }
}

template <class DataTypes, class ElementType, class ExecutionPolicy>
struct ExecPolicyComputeDisplacementStrategy : public ComputeElementForceStrategy<DataTypes, ElementType>
{
    using trait = elasticity::trait<DataTypes, ElementType>;
    using TopologyElement = typename trait::TopologyElement;
    using ElementStiffness = typename trait::ElementStiffness;
    using Real = sofa::Real_t<DataTypes>;

    void compute(
        const sofa::type::vector<TopologyElement>& elements,
        const sofa::VecCoord_t<DataTypes>& position,
        const sofa::VecCoord_t<DataTypes>& restPosition,
        sofa::type::vector<sofa::type::Vec<trait::NumberOfDofsInElement, Real>>& elementForces) final
    {
        if (!m_elementStiffnesses)
        {
            dmsg_error("Compute force") << "Invalid stiffness matrix vector.";
            return;
        }

        if (this->m_elementStiffnesses->size() != elements.size())
        {
            msg_error("Compute force") << "The number of element stiffness matrices is different from the number of elements";
            return;
        }

        std::ranges::iota_view indices {static_cast<decltype(elements.size())>(0ul), elements.size()};

        std::for_each(ExecutionPolicy{}, indices.begin(), indices.end(),
            [&](const auto elementId)
            {
                const auto& element = elements[elementId];
                const auto& stiffnessMatrix = (*this->m_elementStiffnesses)[elementId];

                typename trait::ElementDisplacement displacement{ sofa::type::NOINIT };

                for (sofa::Size j = 0; j < trait::NumberOfNodesInElement; ++j)
                {
                    for (sofa::Size k = 0; k < trait::spatial_dimensions; ++k)
                    {
                        displacement[j * trait::spatial_dimensions + k] = position[element[j]][k] - restPosition[element[j]][k];
                    }
                }

                elementForces[elementId] = stiffnessMatrix * displacement;
            });
    }

    void setElementStiffnessMatrices(const sofa::type::vector<ElementStiffness>& m_elementStiffness) override
    {
        m_elementStiffnesses = &m_elementStiffness;
    }

    sofa::type::vector<ElementStiffness> const* m_elementStiffnesses { nullptr };
};

template <class DataTypes, class ElementType>
void ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::addForce(
    const sofa::core::MechanicalParams* mparams,
    sofa::DataVecDeriv_t<DataTypes>& f,
    const sofa::DataVecCoord_t<DataTypes>& x,
    const sofa::DataVecDeriv_t<DataTypes>& v)
{
    SOFA_UNUSED(mparams);
    SOFA_UNUSED(v);

    auto forceAccessor = sofa::helper::getWriteOnlyAccessor(f);
    auto positionAccessor = sofa::helper::getReadAccessor(x);
    auto restPositionAccessor = this->mstate->readRestPositions();

    if (l_topology == nullptr) return;

    const auto& elements = trait::FiniteElement::getElementSequence(*l_topology);

    m_elementForce.resize(elements.size());

    if (m_computeElementForceStrategy)
    {
        m_computeElementForceStrategy->setElementStiffnessMatrices(this->m_elementStiffness);
        m_computeElementForceStrategy->compute(elements, positionAccessor.ref(), restPositionAccessor.ref(), m_elementForce);
    }

    // dispatch the element force to the degrees of freedom.
    // this operation is done outside the compute strategy because it is not thread-safe.
    for (sofa::Size i = 0; i < elements.size(); ++i)
    {
        const auto& element = elements[i];
        const auto& elementForce = m_elementForce[i];

        for (sofa::Size j = 0; j < trait::NumberOfNodesInElement; ++j)
        {
            auto& nodeForce = forceAccessor[element[j]];
            for (sofa::Size k = 0; k < trait::spatial_dimensions; ++k)
            {
                nodeForce[k] -= elementForce[j * trait::spatial_dimensions + k];
            }
        }
    }
}

template <class DataTypes, class ElementType>
void ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::addDForce(
    const sofa::core::MechanicalParams* mparams,
    sofa::DataVecDeriv_t<DataTypes>& df,
    const sofa::DataVecDeriv_t<DataTypes>& dx)
{
    if (this->isComponentStateInvalid())
        return;

    auto dfAccessor = sofa::helper::getWriteAccessor(df);
    auto dxAccessor = sofa::helper::getReadAccessor(dx);
    dfAccessor.resize(dxAccessor.size());

    const auto& elements = trait::FiniteElement::getElementSequence(*l_topology);

    const auto kFactor = static_cast<sofa::Real_t<DataTypes>>(sofa::core::mechanicalparams::kFactorIncludingRayleighDamping(
            mparams, this->rayleighStiffness.getValue()));

    m_elementDForce.resize(elements.size());

    std::ranges::iota_view indices {static_cast<decltype(elements.size())>(0ul), elements.size()};
    std::for_each(std::execution::par, indices.begin(), indices.end(),
        [&](const auto elementId)
        {
            const auto& element = elements[elementId];
            const auto& stiffnessMatrix = this->m_elementStiffness[elementId];

            sofa::type::Vec<trait::NumberOfDofsInElement, sofa::Real_t<DataTypes>> element_dx(sofa::type::NOINIT);
            for (sofa::Size i = 0; i < trait::NumberOfNodesInElement; ++i)
            {
                VecView<trait::spatial_dimensions, sofa::Real_t<DataTypes>> node_dx(element_dx, i * trait::spatial_dimensions);
                node_dx = dxAccessor[element[i]];
            }

            m_elementDForce[elementId] = (-kFactor) * (stiffnessMatrix * element_dx);
        });

    // dispatch the element dforce to the degrees of freedom.
    // this operation is done outside the compute strategy because it is not thread-safe.
    for (std::size_t elementId = 0; elementId < elements.size(); ++elementId)
    {
        const auto& element = elements[elementId];

        for (sofa::Size i = 0; i < trait::NumberOfNodesInElement; ++i)
        {
            VecView<trait::spatial_dimensions, sofa::Real_t<DataTypes>> nodedForce(m_elementDForce[elementId], i * trait::spatial_dimensions);
            dfAccessor[element[i]] += nodedForce.toVec();
        }
    }
}

template <class DataTypes, class ElementType>
void ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::buildStiffnessMatrix(
    sofa::core::behavior::StiffnessMatrix* matrix)
{
    if (this->isComponentStateInvalid())
        return;

    auto dfdx = matrix->getForceDerivativeIn(this->mstate)
        .withRespectToPositionsIn(this->mstate);

    sofa::type::Mat<trait::spatial_dimensions, trait::spatial_dimensions, sofa::Real_t<DataTypes>> localMatrix(sofa::type::NOINIT);

    const auto& elements = trait::FiniteElement::getElementSequence(*l_topology);
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

    const auto& elements = trait::FiniteElement::getElementSequence(*l_topology);
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

template <class DataTypes, class ElementType>
void ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::selectStrategy()
{
    const std::string& computeStrategy = d_computeStrategy.getValue().getSelectedItem();

    if (computeStrategy == sequencedComputeStrategy)
    {
        m_computeElementForceStrategy = std::make_unique<
            ExecPolicyComputeDisplacementStrategy<DataTypes, ElementType,
        std::execution::sequenced_policy>>();
    }
    else if (computeStrategy == unsequencedComputeStrategy)
    {
        m_computeElementForceStrategy = std::make_unique<
            ExecPolicyComputeDisplacementStrategy<DataTypes, ElementType,
        std::execution::unsequenced_policy>>();
    }
    else if (computeStrategy == parallelComputeStrategy)
    {
        m_computeElementForceStrategy = std::make_unique<
            ExecPolicyComputeDisplacementStrategy<DataTypes, ElementType,
        std::execution::parallel_policy>>();
    }
    else if (computeStrategy == parallelUnsequencedComputeStrategy)
    {
        m_computeElementForceStrategy = std::make_unique<
            ExecPolicyComputeDisplacementStrategy<DataTypes, ElementType,
        std::execution::unsequenced_policy>>();
    }
    else
    {
        msg_error() << "Unknown compute strategy '" + computeStrategy + "'";
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
    }
}

}  // namespace elasticity

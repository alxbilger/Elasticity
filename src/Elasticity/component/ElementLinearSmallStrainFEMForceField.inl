#pragma once
#include <Elasticity/component/ElementLinearSmallStrainFEMForceField.h>
#include <Elasticity/impl/LameParameters.h>
#include <Elasticity/impl/VecView.h>
#include <Elasticity/impl/VectorTools.h>
#include <sofa/core/behavior/ForceField.inl>
#include <ranges>
#include <execution>

namespace elasticity
{

constexpr std::string_view sequentialComputeStrategy = "sequential";
constexpr std::string_view parallelComputeStrategy = "parallel";

template <class DataTypes, class ElementType>
ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::ElementLinearSmallStrainFEMForceField()
    : d_computeStrategy(initData(&d_computeStrategy, "computeStrategy", "The compute strategy used to compute the forces"))
{
    sofa::helper::OptionsGroup computeStrategyOptions{sequentialComputeStrategy, parallelComputeStrategy};
    computeStrategyOptions.setSelectedItem(std::string(parallelComputeStrategy));
    d_computeStrategy.setValue(computeStrategyOptions);

    this->addUpdateCallback("precomputeStiffness", {&this->d_youngModulus, &this->d_poissonRatio},
    [this](const sofa::core::DataTracker& )
    {
        precomputeElementStiffness();
        return this->getComponentState();
    }, {});

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
    sofa::core::behavior::ForceField<DataTypes>::init();

    if (!this->isComponentStateInvalid())
    {
        this->validateTopology();
    }

    if (!this->isComponentStateInvalid())
    {
        selectStrategy();
    }

    if (!this->isComponentStateInvalid())
    {
        precomputeElementStiffness();
    }

    if (!this->isComponentStateInvalid())
    {
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Valid);
    }
}

template <class DataTypes, class ElementType>
struct SequentialComputeDisplacementStrategy : public ComputeElementForceStrategy<DataTypes, ElementType>
{
    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;
    static constexpr sofa::Size NumberOfNodesInElement = ElementType::NumberOfNodes;
    static constexpr sofa::Size NumberOfDofsInElement = NumberOfNodesInElement * spatial_dimensions;
    using TopologyElement = ComputeElementForceStrategy<DataTypes, ElementType>::TopologyElement;
    using ElementDisplacement = ComputeElementForceStrategy<DataTypes, ElementType>::ElementDisplacement;
    using Real = sofa::Real_t<DataTypes>;
    using VecCoord = sofa::VecCoord_t<DataTypes>;
    using ElementStiffness = elasticity::ElementStiffness<DataTypes, ElementType>;

    void compute(
        const sofa::type::vector<TopologyElement>& elements,
        const VecCoord& position,
        const VecCoord& restPosition,
        sofa::type::vector<sofa::type::Vec<NumberOfDofsInElement, Real>>& elementForces) override
    {
        computeElementForceLoop(std::execution::unseq, elements, position, restPosition, elementForces);
    }

    template <class ExecutionPolicy>
    void computeElementForceLoop(
        ExecutionPolicy policy,
        const sofa::type::vector<TopologyElement>& elements,
        const VecCoord& position,
        const VecCoord& restPosition,
        sofa::type::vector<sofa::type::Vec<NumberOfDofsInElement, Real>>& elementForces)
    {
        if (this->m_elementStiffnesses->size() != elements.size())
        {
            msg_error("Compute force") << "The number of element stiffness matrices is different from the number of elements";
            return;
        }

        std::ranges::iota_view indices { static_cast<decltype(elements.size())>(0), elements.size()};

        std::for_each(policy, indices.begin(), indices.end(),
            [&](const auto& id)
            {
                this->computeElementForce(id, elements, position, restPosition, elementForces);
            });
    }

    void computeElementForce(
        sofa::Size elementId,
        const sofa::type::vector<TopologyElement>& elements,
        const VecCoord& position,
        const VecCoord& restPosition,
        sofa::type::vector<sofa::type::Vec<NumberOfDofsInElement, Real>>& elementForces)
    {
        const auto& element = elements[elementId];
        const auto& stiffnessMatrix = (*this->m_elementStiffnesses)[elementId];

        ElementDisplacement displacement{ sofa::type::NOINIT };

        for (sofa::Size j = 0; j < NumberOfNodesInElement; ++j)
        {
            for (sofa::Size k = 0; k < spatial_dimensions; ++k)
            {
                displacement[j * spatial_dimensions + k] = position[element[j]][k] - restPosition[element[j]][k];
            }
        }

        elementForces[elementId] = stiffnessMatrix * displacement;
    }

    void setElementStiffnessMatrices(const sofa::type::vector<ElementStiffness>& m_elementStiffness) override
    {
        m_elementStiffnesses = &m_elementStiffness;
    }

    sofa::type::vector<ElementStiffness> const* m_elementStiffnesses { nullptr };
};

template <class DataTypes, class ElementType>
struct ParallelComputeDisplacementStrategy : public SequentialComputeDisplacementStrategy<DataTypes, ElementType>
{
    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;
    static constexpr sofa::Size NumberOfNodesInElement = ElementType::NumberOfNodes;
    static constexpr sofa::Size NumberOfDofsInElement = NumberOfNodesInElement * spatial_dimensions;
    using TopologyElement = ComputeElementForceStrategy<DataTypes, ElementType>::TopologyElement;
    using Real = sofa::Real_t<DataTypes>;
    using VecCoord = sofa::VecCoord_t<DataTypes>;

    void compute(
        const sofa::type::vector<TopologyElement>& elements,
        const VecCoord& position,
        const VecCoord& restPosition,
        sofa::type::vector<sofa::type::Vec<NumberOfDofsInElement, Real>>& elementForces) override
    {
        this->computeElementForceLoop(std::execution::par_unseq, elements, position, restPosition, elementForces);
    }
};

template <class DataTypes, class ElementType>
void ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::addForce(
    const sofa::core::MechanicalParams* mparams, DataVecDeriv& f, const DataVecCoord& x,
    const DataVecDeriv& v)
{
    SOFA_UNUSED(mparams);
    SOFA_UNUSED(v);

    auto forceAccessor = sofa::helper::getWriteOnlyAccessor(f);
    auto positionAccessor = sofa::helper::getReadAccessor(x);
    auto restPositionAccessor = this->mstate->readRestPositions();

    if (l_topology == nullptr) return;

    const auto& elements = FiniteElement::getElementSequence(*l_topology);

    m_elementForce.resize(elements.size());

    if (m_computeElementForceStrategy)
        m_computeElementForceStrategy->compute(elements, positionAccessor.ref(), restPositionAccessor.ref(), m_elementForce);

    // dispatch the element force to the degrees of freedom
    for (sofa::Size i = 0; i < elements.size(); ++i)
    {
        const auto& element = elements[i];
        const auto& elementForce = m_elementForce[i];

        for (sofa::Size j = 0; j < NumberOfNodesInElement; ++j)
        {
            auto& f_j = forceAccessor[element[j]];
            for (sofa::Size k = 0; k < spatial_dimensions; ++k)
            {
                f_j[k] -= elementForce[j * spatial_dimensions + k];
            }
        }
    }
}

template <class DataTypes, class ElementType>
void ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::addDForce(
    const sofa::core::MechanicalParams* mparams, DataVecDeriv& df, const DataVecDeriv& dx)
{
    if (this->isComponentStateInvalid())
        return;

    auto dfAccessor = sofa::helper::getWriteAccessor(df);
    auto dxAccessor = sofa::helper::getReadAccessor(dx);
    dfAccessor.resize(dxAccessor.size());

    const auto& elements = FiniteElement::getElementSequence(*l_topology);

    const Real kFactor = static_cast<Real>(sofa::core::mechanicalparams::kFactorIncludingRayleighDamping(
            mparams, this->rayleighStiffness.getValue()));

    auto elementStiffnessIt = m_elementStiffness.begin();
    for (const auto& element : elements)
    {
        sofa::type::Vec<NumberOfDofsInElement, Real> element_dx(sofa::type::NOINIT);

        for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
        {
            VecView<spatial_dimensions, Real> node_dx(element_dx, i * spatial_dimensions);
            node_dx = dxAccessor[element[i]];
        }

        const auto& stiffnessMatrix = *elementStiffnessIt++;
        auto dForce = (-kFactor) * (stiffnessMatrix * element_dx);

        for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
        {
            VecView<spatial_dimensions, Real> nodedForce(dForce, i * spatial_dimensions);
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

    sofa::type::Mat<spatial_dimensions, spatial_dimensions, Real> localMatrix(sofa::type::NOINIT);

    const auto& elements = FiniteElement::getElementSequence(*l_topology);
    auto elementStiffnessIt = m_elementStiffness.begin();
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
SReal ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::getPotentialEnergy(
    const sofa::core::MechanicalParams*, const DataVecCoord& x) const
{
    return 0;
}

template <class DataTypes, class ElementType>
void ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::addKToMatrix(
    sofa::linearalgebra::BaseMatrix* matrix, SReal kFact, unsigned& offset)
{
    if (this->isComponentStateInvalid())
        return;

    using LocalMatType = sofa::type::Mat<spatial_dimensions, spatial_dimensions, Real>;
    LocalMatType localMatrix{sofa::type::NOINIT};

    const auto& elements = FiniteElement::getElementSequence(*l_topology);
    auto elementStiffnessIt = m_elementStiffness.begin();
    for (const auto& element : elements)
    {
        const auto& stiffnessMatrix = *elementStiffnessIt++;

        for (sofa::Index n1 = 0; n1 < NumberOfNodesInElement; ++n1)
        {
            for (sofa::Index n2 = 0; n2 < NumberOfNodesInElement; ++n2)
            {
                stiffnessMatrix.getsub(spatial_dimensions * n1, spatial_dimensions * n2, localMatrix); //extract the submatrix corresponding to the coupling of nodes n1 and n2

                const auto value = (-static_cast<Real>(kFact)) * static_cast<ScalarOrMatrix<LocalMatType>>(localMatrix);
                matrix->add(
                   offset + element[n1] * spatial_dimensions,
                   offset + element[n2] * spatial_dimensions, value);
            }
        }
    }
}

template <class DataTypes, class ElementType>
void ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::precomputeElementStiffness()
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

    const auto& elements = FiniteElement::getElementSequence(*l_topology);
    m_elementStiffness.reserve(elements.size());

    for (const auto& element : elements)
    {
        const std::array<Coord, NumberOfNodesInElement> nodesCoordinates = extractNodesVectorFromGlobalVector(element, restPositionAccessor.ref());
        ElementStiffness K = integrate<DataTypes, ElementType>(nodesCoordinates, m_elasticityTensor);
        m_elementStiffness.push_back(K);
    }

    if (m_computeElementForceStrategy)
        m_computeElementForceStrategy->setElementStiffnessMatrices(m_elementStiffness);
    else
        dmsg_error() << "The compute strategy is not yet defined";

    // precompute strain-displacement tensors at the nodes positions
    m_strainDisplacement.clear();
    m_strainDisplacement.resize(elements.size());

    static const std::array<ReferenceCoord, NumberOfNodesInElement>& referenceElementNodes =
        FiniteElement::referenceElementNodes;

    for (sofa::Size i = 0; i < elements.size(); ++i)
    {
        const auto& element = elements[i];
        const auto nodesCoordinates = extractNodesVectorFromGlobalVector(element, restPositionAccessor.ref());
        for (sofa::Size j = 0; j < NumberOfNodesInElement; ++j)
        {
            const ReferenceCoord& x = referenceElementNodes[j];

            // gradient of shape functions in the reference element evaluated at the quadrature point
            const sofa::type::Mat<NumberOfNodesInElement, TopologicalDimension, Real> dN_dq_ref =
                FiniteElement::gradientShapeFunctions(x);

            // jacobian of the mapping from the reference space to the physical space, evaluated at the
            // quadrature point
            sofa::type::Mat<spatial_dimensions, TopologicalDimension, Real> jacobian;
            for (sofa::Size n = 0; n < NumberOfNodesInElement; ++n)
                jacobian += sofa::type::dyad(nodesCoordinates[n], dN_dq_ref[n]);

            const sofa::type::Mat<TopologicalDimension, spatial_dimensions, Real> J_inv =
                elasticity::inverse(jacobian);

            // gradient of the shape functions in the physical element evaluated at the quadrature point
            sofa::type::Mat<NumberOfNodesInElement, spatial_dimensions, Real> dN_dq(sofa::type::NOINIT);
            for (sofa::Size n = 0; n < NumberOfNodesInElement; ++n)
                dN_dq[n] = J_inv.transposed() * dN_dq_ref[n];

            const auto B = makeStrainDisplacement<DataTypes, ElementType>(dN_dq);

            m_strainDisplacement[i][j] = B;
        }
    }
}

template <class DataTypes, class ElementType>
auto ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::computeElementDisplacement(
    const std::array<Coord, NumberOfNodesInElement>& elementNodesCoordinates,
    const std::array<Coord, NumberOfNodesInElement>& restElementNodesCoordinates) -> ElementDisplacement
{
    ElementDisplacement displacement(sofa::type::NOINIT);
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
void ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::selectStrategy()
{
    const std::string& computeStrategy = d_computeStrategy.getValue().getSelectedItem();

    if (computeStrategy == sequentialComputeStrategy)
    {
        m_computeElementForceStrategy = std::make_unique<SequentialComputeDisplacementStrategy<DataTypes, ElementType>>();
    }
    else if (computeStrategy == parallelComputeStrategy)
    {
        m_computeElementForceStrategy = std::make_unique<ParallelComputeDisplacementStrategy<DataTypes, ElementType>>();
    }
    else
    {
        msg_error() << "Unknown compute strategy '" + computeStrategy + "'";
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
    }
}

}  // namespace elasticity

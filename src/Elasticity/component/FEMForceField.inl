#pragma once
#include <Elasticity/component/FEMForceField.h>
#include <Elasticity/impl/VecView.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/simulation/task/ParallelForEach.h>

namespace elasticity
{

template <class DataTypes, class ElementType>
FEMForceField<DataTypes, ElementType>::FEMForceField()
    : d_computeForceStrategy(initData(&d_computeForceStrategy, "computeForceStrategy", std::string("The compute strategy used to compute the forces.\n" + ComputeStrategy::dataDescription()).c_str()))
    , d_computeForceDerivStrategy(initData(&d_computeForceDerivStrategy, "computeForceDerivStrategy", std::string("The compute strategy used to compute the forces derivatives.\n" + ComputeStrategy::dataDescription()).c_str()))
    , d_elementSpace(initData(&d_elementSpace, static_cast<sofa::Real_t<DataTypes>>(0.125), "elementSpace", "When rendering, the space between elements"))
{
    d_elementSpace.setGroup("Visualization");
}

template <class DataTypes, class ElementType>
void FEMForceField<DataTypes, ElementType>::init()
{
    sofa::core::behavior::ForceField<DataTypes>::init();

    if (!this->isComponentStateInvalid())
    {
        TopologyAccessor::init();
    }

    if (!this->isComponentStateInvalid())
    {
        this->initTaskScheduler();
    }
}

template <class DataTypes, class ElementType>
void FEMForceField<DataTypes, ElementType>::addForce(
    const sofa::core::MechanicalParams* mparams,
    sofa::DataVecDeriv_t<DataTypes>& f,
    const sofa::DataVecCoord_t<DataTypes>& x,
    const sofa::DataVecDeriv_t<DataTypes>& v)
{
    SOFA_UNUSED(mparams);
    SOFA_UNUSED(v);

    auto positionAccessor = sofa::helper::getReadAccessor(x);

    if (this->l_topology == nullptr) return;

    const auto& elements = trait::FiniteElement::getElementSequence(*this->l_topology);
    m_elementForce.resize(elements.size());

    this->computeElementsForces(mparams, m_elementForce, positionAccessor.ref());

    auto forceAccessor = sofa::helper::getWriteOnlyAccessor(f);

    // dispatch the element force to the degrees of freedom.
    // this operation is done outside the compute strategy because it is not thread-safe.
    dispatchElementForcesToNodes(elements, forceAccessor.wref());
}

template <class DataTypes, class ElementType>
void FEMForceField<DataTypes, ElementType>::computeElementsForces(
    const sofa::core::MechanicalParams* mparams,
    sofa::type::vector<ElementForce>& f,
    const sofa::VecCoord_t<DataTypes>& x)
{
    this->beforeElementForce(mparams, f, x);

    const auto& elements = trait::FiniteElement::getElementSequence(*this->l_topology);

    auto computeForceStrategyAccessor = sofa::helper::getReadAccessor(d_computeForceStrategy);
    const auto& computeForceStrategy = computeForceStrategyAccessor->key();

    sofa::simulation::ForEachExecutionPolicy executionPolicy(sofa::simulation::ForEachExecutionPolicy::SEQUENTIAL);
    if (computeForceStrategy == parallelComputeStrategy)
        executionPolicy = sofa::simulation::ForEachExecutionPolicy::PARALLEL;

    sofa::simulation::forEachRange(executionPolicy, *this->m_taskScheduler,
        static_cast<decltype(elements.size())>(0), elements.size(), [this, mparams, &f, &x](const auto& range)
        {
            this->computeElementsForces(range, mparams, f, x);
        });
}

template <class DataTypes, class ElementType>
void FEMForceField<DataTypes, ElementType>::dispatchElementForcesToNodes(
    const sofa::type::vector<typename trait::TopologyElement>& elements,
    sofa::VecDeriv_t<DataTypes>& nodeForces)
{
    for (sofa::Size i = 0; i < elements.size(); ++i)
    {
        const auto& element = elements[i];
        const auto& elementForce = m_elementForce[i];

        for (sofa::Size j = 0; j < trait::NumberOfNodesInElement; ++j)
        {
            auto& nodeForce = nodeForces[element[j]];
            for (sofa::Size k = 0; k < trait::spatial_dimensions; ++k)
            {
                nodeForce[k] -= elementForce[j * trait::spatial_dimensions + k];
            }
        }
    }
}


template <class DataTypes, class ElementType>
void FEMForceField<DataTypes, ElementType>::addDForce(
    const sofa::core::MechanicalParams* mparams,
    sofa::DataVecDeriv_t<DataTypes>& df,
    const sofa::DataVecDeriv_t<DataTypes>& dx)
{
    if (this->isComponentStateInvalid())
        return;

    auto dfAccessor = sofa::helper::getWriteAccessor(df);
    auto dxAccessor = sofa::helper::getReadAccessor(dx);
    dfAccessor.resize(dxAccessor.size());

    const auto& elements = trait::FiniteElement::getElementSequence(*this->l_topology);

    const auto kFactor = static_cast<sofa::Real_t<DataTypes>>(sofa::core::mechanicalparams::kFactorIncludingRayleighDamping(
            mparams, this->rayleighStiffness.getValue()));

    m_elementDForce.resize(elements.size());

    this->computeElementsForcesDeriv(mparams, m_elementDForce, dxAccessor.ref(), kFactor);

    // dispatch the element dforce to the degrees of freedom.
    // this operation is done outside the compute strategy because it is not thread-safe.
    for (std::size_t elementId = 0; elementId < elements.size(); ++elementId)
    {
        const auto& element = elements[elementId];

        for (sofa::Size i = 0; i < trait::NumberOfNodesInElement; ++i)
        {
            VecView<trait::spatial_dimensions, sofa::Real_t<DataTypes>> nodedForce(m_elementDForce[elementId], i * trait::spatial_dimensions);
            dfAccessor[element[i]] -= nodedForce.toVec();
        }
    }
}

template <class DataTypes, class ElementType>
void FEMForceField<DataTypes, ElementType>::computeElementsForcesDeriv(
    const sofa::core::MechanicalParams* mparams, sofa::type::vector<ElementForce>& df,
    const sofa::VecDeriv_t<DataTypes>& dx, sofa::Real_t<DataTypes> kFactor)
{
    const auto& elements = trait::FiniteElement::getElementSequence(*this->l_topology);

    auto computeForceDerivStrategyAccessor = sofa::helper::getReadAccessor(d_computeForceDerivStrategy);
    const auto& computeForceDerivStrategy = computeForceDerivStrategyAccessor->key();

    sofa::simulation::ForEachExecutionPolicy executionPolicy(sofa::simulation::ForEachExecutionPolicy::SEQUENTIAL);
    if (computeForceDerivStrategy == parallelComputeStrategy)
        executionPolicy = sofa::simulation::ForEachExecutionPolicy::PARALLEL;

    sofa::simulation::forEachRange(executionPolicy, *this->m_taskScheduler,
        static_cast<decltype(elements.size())>(0), elements.size(), [this, mparams, &df, &dx, kFactor](const auto& range)
        {
            this->computeElementsForcesDeriv(range, mparams, df, dx, kFactor);
        });
}

template <class DataTypes, class ElementType>
void FEMForceField<DataTypes, ElementType>::draw(const sofa::core::visual::VisualParams* vparams)
{
    if (!vparams->displayFlags().getShowForceFields())
        return;

    if (!this->l_topology)
        return;

    const auto stateLifeCycle = vparams->drawTool()->makeStateLifeCycle();

    if (vparams->displayFlags().getShowWireFrame())
        vparams->drawTool()->setPolygonMode(0, true);

    const auto& x = this->mstate->read(sofa::core::vec_id::read_access::position)->getValue();

    if constexpr (std::same_as<ElementType, sofa::geometry::Triangle> || std::same_as<ElementType, sofa::geometry::Tetrahedron> || std::same_as<ElementType, sofa::geometry::Hexahedron> )
    {
        m_drawMesh.elementSpace = d_elementSpace.getValue();
        m_drawMesh.drawAllElements(vparams->drawTool(), x, this->l_topology.get());
    }
}

}  // namespace elasticity

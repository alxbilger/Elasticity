#pragma once
#include <Elasticity/component/LinearSmallStrainFEMForceField.h>
#include <Elasticity/finiteelement/FiniteElement[Edge].h>
#include <Elasticity/finiteelement/FiniteElement[Hexahedron].h>
#include <Elasticity/finiteelement/FiniteElement[Quad].h>
#include <Elasticity/finiteelement/FiniteElement[Tetrahedron].h>
#include <Elasticity/finiteelement/FiniteElement[Triangle].h>

#include <sofa/core/ObjectFactory.h>

namespace elasticity
{

template <class DataTypes>
LinearSmallStrainFEMForceField<DataTypes>::LinearSmallStrainFEMForceField()
    : l_topology(initLink("topology", "Link to a topology containing elements"))
    , d_poissonRatio(initData(&d_poissonRatio, static_cast<Real>(0.45), "poissonRatio", "Poisson's ratio"))
    , d_youngModulus(initData(&d_youngModulus, static_cast<Real>(1e6), "youngModulus", "Young's modulus"))
    , d_vonMisesStressValues(initData(&d_vonMisesStressValues, "vonMisesStressValues", "Von Mises stress values"))
{
    static std::string groupName = "Mechanical parameter";
    d_poissonRatio.setGroup(groupName);
    d_youngModulus.setGroup(groupName);

    d_vonMisesStressValues.setGroup("Output");
}

template <class DataTypes>
void LinearSmallStrainFEMForceField<DataTypes>::init()
{
    sofa::core::behavior::ForceField<DataTypes>::init();

    if (!this->isComponentStateInvalid())
    {
        this->validateTopology();
    }

    if (!this->isComponentStateInvalid())
    {
        selectFEMTypes();
    }

    if (!this->isComponentStateInvalid())
    {
        precomputeElementStiffness();
    }

    if (!this->isComponentStateInvalid())
    {
        m_vonMisesStressContainer.resize(this->mstate->getSize());
        d_vonMisesStressValues.setValue(m_vonMisesStressContainer.getStressValues());
    }

    if (!this->isComponentStateInvalid())
    {
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Valid);
    }
}

template <class DataTypes>
void LinearSmallStrainFEMForceField<DataTypes>::addForce(
    const sofa::core::MechanicalParams* mparams, DataVecDeriv& f, const DataVecCoord& x,
    const DataVecDeriv& v)
{
    SOFA_UNUSED(mparams);
    SOFA_UNUSED(v);

    auto forceAccessor = sofa::helper::getWriteOnlyAccessor(f);
    auto positionAccessor = sofa::helper::getReadAccessor(x);
    auto restPositionAccessor = this->mstate->readRestPositions();

    m_vonMisesStressContainer.clear();
    m_vonMisesStressContainer.resize(positionAccessor.size());

    for (const auto& finiteElement : m_finiteElements)
    {
        finiteElement->addForce(forceAccessor.wref(), positionAccessor.ref(), restPositionAccessor.ref());
    }

    for (const auto& finiteElement : m_finiteElements)
    {
        finiteElement->computeVonMisesStress(m_vonMisesStressContainer, positionAccessor.ref(), restPositionAccessor.ref());
    }

    m_vonMisesStressContainer.accept();

    d_vonMisesStressValues.setValue(m_vonMisesStressContainer.getStressValues());
}

template <class DataTypes>
void LinearSmallStrainFEMForceField<DataTypes>::addDForce(
    const sofa::core::MechanicalParams* mparams, DataVecDeriv& df, const DataVecDeriv& dx)
{
    auto dfAccessor = sofa::helper::getWriteAccessor(df);
    auto dxAccessor = sofa::helper::getReadAccessor(dx);
    dfAccessor.resize(dxAccessor.size());

    const Real kFactor = (Real)sofa::core::mechanicalparams::kFactorIncludingRayleighDamping(
        mparams, this->rayleighStiffness.getValue());

    for (const auto& finiteElement : m_finiteElements)
    {
        finiteElement->addDForce(dfAccessor.wref(), dxAccessor.ref(), kFactor);
    }
}

template <class DataTypes>
void LinearSmallStrainFEMForceField<DataTypes>::buildStiffnessMatrix(
    sofa::core::behavior::StiffnessMatrix* matrix)
{
    auto dfdx = matrix->getForceDerivativeIn(this->mstate)
               .withRespectToPositionsIn(this->mstate);

    for (const auto& finiteElement : m_finiteElements)
    {
        finiteElement->buildStiffnessMatrix(dfdx);
    }
}

template <class DataTypes>
SReal LinearSmallStrainFEMForceField<DataTypes>::getPotentialEnergy(
    const sofa::core::MechanicalParams*, const DataVecCoord& x) const
{
    return 0;
}

template <class DataTypes>
void LinearSmallStrainFEMForceField<DataTypes>::validateTopology()
{
    if (l_topology.empty())
    {
        msg_info() << "Link to Topology container should be set to ensure right behavior. First "
                      "Topology found in current context will be used.";
        l_topology.set(this->getContext()->getMeshTopologyLink());
    }

    if (l_topology == nullptr)
    {
        msg_error() << "No topology component found at path: " << this->l_topology.getLinkedPath()
                    << ", nor in current context: " << this->getContext()->name
                    << ". Object must have a BaseMeshTopology. "
                    << "The list of available BaseMeshTopology components is: "
                    << sofa::core::ObjectFactory::getInstance()
                           ->listClassesDerivedFrom<sofa::core::topology::BaseMeshTopology>();
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
    }
}

template <class DataTypes>
void LinearSmallStrainFEMForceField<DataTypes>::precomputeElementStiffness()
{
    auto restPositionAccessor = this->mstate->readRestPositions();
    for (const auto& finiteElement : m_finiteElements)
    {
        finiteElement->precomputeElementStiffness(restPositionAccessor.ref(), d_youngModulus.getValue(), d_poissonRatio.getValue());
    }
}

template <class DataTypes>
void LinearSmallStrainFEMForceField<DataTypes>::selectFEMTypes()
{
    if constexpr (spatial_dimensions == 1)
    {
        addLinearFEMType<sofa::geometry::Edge>();
    }
    else if constexpr (spatial_dimensions == 2)
    {
        const auto nbTriangles = this->l_topology->getNbTriangles();
        const auto nbQuads = this->l_topology->getNbQuads();
        const auto nbEdges = this->l_topology->getNbEdges();

        if (nbTriangles > 0 || nbQuads > 0)
        {
            addLinearFEMType<sofa::geometry::Triangle>();
            addLinearFEMType<sofa::geometry::Quad>();
        }
        else if (nbEdges > 0)
        {
            addLinearFEMType<sofa::geometry::Edge>();
        }
        else
        {
            addLinearFEMType<sofa::geometry::Triangle>();
            addLinearFEMType<sofa::geometry::Quad>();
        }
    }
    else if constexpr (spatial_dimensions == 3)
    {
        const auto nbTetras = this->l_topology->getNbTetrahedra();
        const auto nbHexas = this->l_topology->getNbHexahedra();
        const auto nbTriangles = this->l_topology->getNbTriangles();
        const auto nbQuads = this->l_topology->getNbQuads();
        const auto nbEdges = this->l_topology->getNbEdges();

        if (nbTetras > 0 || nbHexas > 0)
        {
            addLinearFEMType<sofa::geometry::Tetrahedron>();
            addLinearFEMType<sofa::geometry::Hexahedron>();
        }
        else if (nbTriangles > 0 || nbQuads > 0)
        {
            addLinearFEMType<sofa::geometry::Triangle>();
            addLinearFEMType<sofa::geometry::Quad>();
        }
        else if (nbEdges > 0)
        {
            addLinearFEMType<sofa::geometry::Edge>();
        }
        else
        {
            addLinearFEMType<sofa::geometry::Tetrahedron>();
            addLinearFEMType<sofa::geometry::Hexahedron>();
        }
    }
}

template <class DataTypes>
template <class ElementType>
void LinearSmallStrainFEMForceField<DataTypes>::addLinearFEMType()
{
    m_finiteElements.emplace_back(
        std::make_unique<LinearFEM<DataTypes, ElementType>>(this->l_topology.get()));
}

}  // namespace elasticity

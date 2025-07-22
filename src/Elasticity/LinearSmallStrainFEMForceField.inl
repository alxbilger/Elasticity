#pragma once
#include <Elasticity/LinearSmallStrainFEMForceField.h>
#include <Elasticity/BaseLinearSmallStrainFEMForceField.inl>

#include <Elasticity/FiniteElement[Edge].h>
#include <Elasticity/FiniteElement[Hexahedron].h>
#include <Elasticity/FiniteElement[Quad].h>
#include <Elasticity/FiniteElement[Tetrahedron].h>
#include <Elasticity/FiniteElement[Triangle].h>

namespace elasticity
{

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

    for (const auto& finiteElement : m_finiteElements)
    {
        finiteElement->addForce(forceAccessor.wref(), positionAccessor.ref(), restPositionAccessor.ref());
    }
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
        addFEMType<sofa::geometry::Edge>();
    }
    else if constexpr (spatial_dimensions == 2)
    {
        const auto nbTriangles = this->l_topology->getNbTriangles();
        const auto nbQuads = this->l_topology->getNbQuads();
        const auto nbEdges = this->l_topology->getNbEdges();

        if (nbTriangles > 0 || nbQuads > 0)
        {
            addFEMType<sofa::geometry::Triangle>();
            addFEMType<sofa::geometry::Quad>();
        }
        else if (nbEdges > 0)
        {
            addFEMType<sofa::geometry::Edge>();
        }
        else
        {
            addFEMType<sofa::geometry::Triangle>();
            addFEMType<sofa::geometry::Quad>();
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
            addFEMType<sofa::geometry::Tetrahedron>();
            addFEMType<sofa::geometry::Hexahedron>();
        }
        else if (nbTriangles > 0 || nbQuads > 0)
        {
            addFEMType<sofa::geometry::Triangle>();
            addFEMType<sofa::geometry::Quad>();
        }
        else if (nbEdges > 0)
        {
            addFEMType<sofa::geometry::Edge>();
        }
        else
        {
            addFEMType<sofa::geometry::Tetrahedron>();
            addFEMType<sofa::geometry::Hexahedron>();
        }
    }
}

template <class DataTypes>
template <class ElementType>
void LinearSmallStrainFEMForceField<DataTypes>::addFEMType()
{
    m_finiteElements.emplace_back(
        std::make_unique<LinearFEM<DataTypes, ElementType>>(this->l_topology.get()));
}

}  // namespace elasticity

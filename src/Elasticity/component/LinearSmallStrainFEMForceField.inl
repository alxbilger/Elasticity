#pragma once
#include <Elasticity/component/LinearSmallStrainFEMForceField.h>
#include <Elasticity/component/BaseElasticityFEMForceField.inl>
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
    : d_poissonRatio(initData(&d_poissonRatio, static_cast<Real>(0.45), "poissonRatio", "Poisson's ratio"))
    , d_youngModulus(initData(&d_youngModulus, static_cast<Real>(1e6), "youngModulus", "Young's modulus"))
    , d_computeVonMisesStress(initData(&d_computeVonMisesStress, true, "computeVonMisesStress", "Compute Von Mises stress"))
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
    BaseElasticityFEMForceField<DataTypes>::init();

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
void LinearSmallStrainFEMForceField<DataTypes>::computeVonMisesStress(const DataVecCoord& x)
{
    auto positionAccessor = sofa::helper::getReadAccessor(x);
    auto restPositionAccessor = this->mstate->readRestPositions();

    if (d_computeVonMisesStress.getValue())
    {
        m_vonMisesStressContainer.clear();
        m_vonMisesStressContainer.resize(positionAccessor.size());

        for (const auto& finiteElement : m_finiteElements)
        {
            finiteElement->computeVonMisesStress(m_vonMisesStressContainer, positionAccessor.ref(),
                                                 restPositionAccessor.ref());
        }

        m_vonMisesStressContainer.accept();

        d_vonMisesStressValues.setValue(m_vonMisesStressContainer.getStressValues());
    }
}

template <class DataTypes>
void LinearSmallStrainFEMForceField<DataTypes>::addForce(
    const sofa::core::MechanicalParams* mparams, DataVecDeriv& f, const DataVecCoord& x,
    const DataVecDeriv& v)
{
    BaseElasticityFEMForceField<DataTypes>::addForce(mparams, f, x, v);
    computeVonMisesStress(x);
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

template <class DataTypes>
void LinearSmallStrainFEMForceField<DataTypes>::applyLambda(
    const std::function<void(BaseFEM<DataTypes>&)>& callable)
{
    for (const auto& finiteElement : m_finiteElements)
    {
        callable(*finiteElement);
    }
}

}  // namespace elasticity

#pragma once

#include <Elasticity/component/HyperelasticityFEMForceField.h>
#include <Elasticity/component/BaseElasticityFEMForceField.inl>
#include <Elasticity/impl/FEM/NonLinearFEM.inl>

namespace elasticity
{

template <class DataTypes>
SReal HyperelasticityFEMForceField<DataTypes>::getPotentialEnergy(
    const sofa::core::MechanicalParams*, const DataVecCoord& x) const
{
    return 0;
}

template <class DataTypes>
void HyperelasticityFEMForceField<DataTypes>::selectFEMTypes()
{
    if constexpr (spatial_dimensions == 1)
    {
        addNonLinearFEMType<sofa::geometry::Edge>();
    }
    else if constexpr (spatial_dimensions == 2)
    {
        const auto nbTriangles = this->l_topology->getNbTriangles();
        const auto nbQuads = this->l_topology->getNbQuads();
        const auto nbEdges = this->l_topology->getNbEdges();

        if (nbTriangles > 0 || nbQuads > 0)
        {
            addNonLinearFEMType<sofa::geometry::Triangle>();
            addNonLinearFEMType<sofa::geometry::Quad>();
        }
        else if (nbEdges > 0)
        {
            addNonLinearFEMType<sofa::geometry::Edge>();
        }
        else
        {
            addNonLinearFEMType<sofa::geometry::Triangle>();
            addNonLinearFEMType<sofa::geometry::Quad>();
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
            addNonLinearFEMType<sofa::geometry::Tetrahedron>();
            addNonLinearFEMType<sofa::geometry::Hexahedron>();
        }
        else if (nbTriangles > 0 || nbQuads > 0)
        {
            addNonLinearFEMType<sofa::geometry::Triangle>();
            addNonLinearFEMType<sofa::geometry::Quad>();
        }
        else if (nbEdges > 0)
        {
            addNonLinearFEMType<sofa::geometry::Edge>();
        }
        else
        {
            addNonLinearFEMType<sofa::geometry::Tetrahedron>();
            addNonLinearFEMType<sofa::geometry::Hexahedron>();
        }
    }
}

template <class DataTypes>
template <class ElementType>
void HyperelasticityFEMForceField<DataTypes>::addNonLinearFEMType()
{
    m_finiteElements.emplace_back(
        std::make_unique<NonLinearFEM<DataTypes, ElementType>>(this->l_topology.get()));
}

template <class DataTypes>
void HyperelasticityFEMForceField<DataTypes>::applyLambda(
    const std::function<void(BaseFEM<DataTypes>&)>& callable)
{
    for (const auto& finiteElement : m_finiteElements)
    {
        callable(*finiteElement);
    }
}

}  // namespace elasticity

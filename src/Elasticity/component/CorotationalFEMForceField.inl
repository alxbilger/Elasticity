#pragma once

#include <Elasticity/component/CorotationalFEMForceField.h>
#include <Elasticity/component/LinearSmallStrainFEMForceField.inl>
#include <Elasticity/impl/CorotationalFEM.h>

namespace elasticity
{

template <class DataTypes>
void CorotationalFEMForceField<DataTypes>::selectFEMTypes()
{
    if constexpr (spatial_dimensions == 2)
    {
        addCorotationalFEMType<sofa::geometry::Triangle>();
        addCorotationalFEMType<sofa::geometry::Quad>();
    }
    else if constexpr (spatial_dimensions == 3)
    {
        const auto nbTriangles = this->l_topology->getNbTriangles();
        const auto nbQuads = this->l_topology->getNbQuads();
        const auto nbTetras = this->l_topology->getNbTetrahedra();
        const auto nbHexas = this->l_topology->getNbHexahedra();

        if (nbTetras > 0 || nbHexas > 0)
        {
            addCorotationalFEMType<sofa::geometry::Tetrahedron>();
            addCorotationalFEMType<sofa::geometry::Hexahedron>();
        }
        else if (nbTriangles > 0 || nbQuads > 0)
        {
            addCorotationalFEMType<sofa::geometry::Triangle>();
            addCorotationalFEMType<sofa::geometry::Quad>();
        }
        else
        {
            addCorotationalFEMType<sofa::geometry::Tetrahedron>();
            addCorotationalFEMType<sofa::geometry::Hexahedron>();
        }
    }
}

template <class DataTypes>
template <class ElementType>
void CorotationalFEMForceField<DataTypes>::addCorotationalFEMType()
{
    this->m_finiteElements.emplace_back(
        std::make_unique<CorotationalFEM<DataTypes, ElementType>>(this->l_topology.get()));
}

}  // namespace elasticity

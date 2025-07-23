#pragma once

#include <Elasticity/CorotationalFEMForceField.h>
#include <Elasticity/LinearSmallStrainFEMForceField.inl>
#include <Elasticity/CorotationalFEM.h>

namespace elasticity
{

template <class DataTypes>
void CorotationalFEMForceField<DataTypes>::selectFEMTypes()
{
    this->m_finiteElements.emplace_back(
        std::make_unique<CorotationalFEM<DataTypes, sofa::geometry::Tetrahedron>>(this->l_topology.get()));
    this->m_finiteElements.emplace_back(
        std::make_unique<CorotationalFEM<DataTypes, sofa::geometry::Hexahedron>>(this->l_topology.get()));
}

}  // namespace elasticity

#pragma once

#include <Elasticity/component/ElementPrefab.h>

namespace elasticity
{

template <class ElementClass>
void ElementPrefab<ElementClass>::init()
{
    TopologyAccessor::init();

    if (!this->isComponentStateInvalid())
    {
        instantiateComponentBasedOnElementTypes();
    }

    if (!this->isComponentStateInvalid())
    {
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Valid);
    }
}

template <class ElementClass>
void ElementPrefab<ElementClass>::instantiateComponentBasedOnElementTypes()
{
    using EdgeType = typename ElementClass::EdgeType;
    using TriangleType = typename ElementClass::TriangleType;
    using QuadType = typename ElementClass::QuadType;
    using TetrahedronType = typename ElementClass::TetrahedronType;
    using HexahedronType = typename ElementClass::HexahedronType;

    if constexpr (spatial_dimensions == 1)
    {
        addComponent<typename ElementClass::EdgeType>();
    }
    else if constexpr (spatial_dimensions == 2)
    {
        const auto nbTriangles = this->l_topology->getNbTriangles();
        const auto nbQuads = this->l_topology->getNbQuads();
        const auto nbEdges = this->l_topology->getNbEdges();

        if (nbTriangles > 0 || nbQuads > 0)
        {
            addComponent<TriangleType>();
            addComponent<QuadType>();
        }
        else if (nbEdges > 0)
        {
            addComponent<EdgeType>();
        }
        else
        {
            msg_warning() << "Cannot find any element in the topology: triangles and quads will be considered by default.";
            addComponent<TriangleType>();
            addComponent<QuadType>();
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
            addComponent<TetrahedronType>();
            addComponent<HexahedronType>();
        }
        else if (nbTriangles > 0 || nbQuads > 0)
        {
            addComponent<TriangleType>();
            addComponent<QuadType>();
        }
        else if (nbEdges > 0)
        {
            addComponent<EdgeType>();
        }
        else
        {
            msg_warning() << "Cannot find any element in the topology: tetrahedron and hexahedron will be considered by default.";
            addComponent<TetrahedronType>();
            addComponent<HexahedronType>();
        }
    }
}

template <class ElementClass>
template <class ComponentType>
void ElementPrefab<ElementClass>::addComponent()
{
    BaseObject::SPtr component = sofa::core::objectmodel::New<ComponentType>();
    assert(component);

    this->getContext()->addObject(component);

    // this will link all the Data with the same name
    component->setSrc("@"+this->getName(), this);

    component->setName(this->getContext()->getNameHelper().resolveName(component->getClassName(), sofa::core::ComponentNameHelper::Convention::python));
}

}  // namespace elasticity

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
            addComponent<typename ElementClass::TriangleType>();
            addComponent<typename ElementClass::QuadType>();
        }
        else if (nbEdges > 0)
        {
            addComponent<typename ElementClass::EdgeType>();
        }
        else
        {
            addComponent<typename ElementClass::TriangleType>();
            addComponent<typename ElementClass::QuadType>();
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
            addComponent<typename ElementClass::TetrahedronType>();
            addComponent<typename ElementClass::HexahedronType>();
        }
        else if (nbTriangles > 0 || nbQuads > 0)
        {
            addComponent<typename ElementClass::TriangleType>();
            addComponent<typename ElementClass::QuadType>();
        }
        else if (nbEdges > 0)
        {
            addComponent<typename ElementClass::EdgeType>();
        }
        else
        {
            addComponent<typename ElementClass::TetrahedronType>();
            addComponent<typename ElementClass::HexahedronType>();
        }
    }
}

template <class ElementClass>
template <class ComponentType>
void ElementPrefab<ElementClass>::addComponent()
{
    BaseObject::SPtr component = sofa::core::objectmodel::New<ComponentType>();
    this->getContext()->addObject(component);

    component->setSrc("@"+this->getName(), this);

    component->setName(this->getContext()->getNameHelper().resolveName(component->getClassName(), sofa::core::ComponentNameHelper::Convention::python));
}

}  // namespace elasticity

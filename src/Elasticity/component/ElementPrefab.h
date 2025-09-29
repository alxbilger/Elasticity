#pragma once

#include <Elasticity/config.h>
#include <sofa/core/objectmodel/BaseObject.h>

#include <Elasticity/component/TopologyAccessor.h>

namespace elasticity
{

/**
 * An intermediate class that instantiates other components based on the types of elements found in
 * the linked topology.
 *
 * @tparam ElementTypes A template providing the types of components to instantiate depending on
 * the type of elements.
 */
template <class ElementTypes>
class ElementPrefab : public virtual TopologyAccessor
{
public:
    SOFA_CLASS(ElementPrefab<ElementTypes>, TopologyAccessor);

private:
    using DataTypes = typename ElementTypes::DataTypes;
    using Real = sofa::Real_t<DataTypes>;
    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;

public:

    void init() override;

protected:

    void instantiateComponentBasedOnElementTypes();

     /**
      * Instantiate a component of type ComponentType and add it to the current context. It also
      * creates a link between Data from this component to the created one.
      */
    template <class ComponentType>
    void addComponent();
};

template<template<class, class> class ElementComponent, class TDataTypes>
struct ElementPrefabTrait
{
    using DataTypes = TDataTypes;

    using EdgeType = ElementComponent<DataTypes, sofa::geometry::Edge>;
    using TriangleType = ElementComponent<DataTypes, sofa::geometry::Triangle>;
    using QuadType = ElementComponent<DataTypes, sofa::geometry::Quad>;
    using TetrahedronType = ElementComponent<DataTypes, sofa::geometry::Tetrahedron>;
    using HexahedronType = ElementComponent<DataTypes, sofa::geometry::Hexahedron>;
};

}  // namespace elasticity

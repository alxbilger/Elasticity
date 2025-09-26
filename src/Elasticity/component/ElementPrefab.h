#pragma once

#include <Elasticity/config.h>
#include <sofa/core/objectmodel/BaseObject.h>

#include <Elasticity/component/TopologyAccessor.h>

namespace elasticity
{

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

    template<class ComponentType>
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

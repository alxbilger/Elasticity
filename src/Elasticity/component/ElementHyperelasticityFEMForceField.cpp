#define ELASTICITY_COMPONENT_ELEMENT_HYPERLASTICITY_FEM_FORCE_FIELD_CPP

#include <Elasticity/component/ElementHyperelasticityFEMForceField.inl>

#include <Elasticity/finiteelement/FiniteElement[all].h>
#include <sofa/core/ObjectFactory.h>

namespace elasticity
{

void registerElementHyperelasticityFEMForceField(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData("Hyperelasticity on linear beams")
        .add< ElementHyperelasticityFEMForceField<sofa::defaulttype::Vec1Types, sofa::geometry::Edge> >()
        .add< ElementHyperelasticityFEMForceField<sofa::defaulttype::Vec2Types, sofa::geometry::Edge> >()
        .add< ElementHyperelasticityFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Edge> >(true));

    factory->registerObjects(sofa::core::ObjectRegistrationData("Hyperelasticity on linear triangles")
        .add< ElementHyperelasticityFEMForceField<sofa::defaulttype::Vec2Types, sofa::geometry::Triangle> >()
        .add< ElementHyperelasticityFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Triangle> >(true));

    factory->registerObjects(sofa::core::ObjectRegistrationData("Hyperelasticity on linear quads")
        .add< ElementHyperelasticityFEMForceField<sofa::defaulttype::Vec2Types, sofa::geometry::Quad> >()
        .add< ElementHyperelasticityFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Quad> >(true));

    factory->registerObjects(sofa::core::ObjectRegistrationData("Hyperelasticity on linear tetrahedra")
        .add< ElementHyperelasticityFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Tetrahedron> >(true));

    factory->registerObjects(sofa::core::ObjectRegistrationData("Hyperelasticity on linear hexahedra")
        .add< ElementHyperelasticityFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Hexahedron> >(true));
}

template class ELASTICITY_API ElementHyperelasticityFEMForceField<sofa::defaulttype::Vec1Types, sofa::geometry::Edge>;
template class ELASTICITY_API ElementHyperelasticityFEMForceField<sofa::defaulttype::Vec2Types, sofa::geometry::Edge>;
template class ELASTICITY_API ElementHyperelasticityFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Edge>;
template class ELASTICITY_API ElementHyperelasticityFEMForceField<sofa::defaulttype::Vec2Types, sofa::geometry::Triangle>;
template class ELASTICITY_API ElementHyperelasticityFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Triangle>;
template class ELASTICITY_API ElementHyperelasticityFEMForceField<sofa::defaulttype::Vec2Types, sofa::geometry::Quad>;
template class ELASTICITY_API ElementHyperelasticityFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Quad>;
template class ELASTICITY_API ElementHyperelasticityFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Tetrahedron>;
template class ELASTICITY_API ElementHyperelasticityFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Hexahedron>;

}

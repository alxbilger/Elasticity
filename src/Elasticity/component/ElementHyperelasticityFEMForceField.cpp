#define ELASTICITY_COMPONENT_ELEMENT_HYPERLASTICITY_FEM_FORCE_FIELD_CPP

#include <Elasticity/component/ElementHyperelasticityFEMForceField.inl>

#include <sofa/fem/FiniteElement[all].h>
#include <sofa/core/ObjectFactory.h>

namespace elasticity
{

void registerElementHyperelasticityFEMForceField(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData("Hyperelasticity")
        .add< ElementHyperelasticityFEMForceField<sofa::defaulttype::Vec1Types, sofa::geometry::Edge> >()
        .add< ElementHyperelasticityFEMForceField<sofa::defaulttype::Vec2Types, sofa::geometry::Edge> >()
        .add< ElementHyperelasticityFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Edge> >()
        .add< ElementHyperelasticityFEMForceField<sofa::defaulttype::Vec2Types, sofa::geometry::Triangle> >()
        .add< ElementHyperelasticityFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Triangle> >()
        .add< ElementHyperelasticityFEMForceField<sofa::defaulttype::Vec2Types, sofa::geometry::Quad> >()
        .add< ElementHyperelasticityFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Quad> >()
        .add< ElementHyperelasticityFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Tetrahedron> >()
        .add< ElementHyperelasticityFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Hexahedron> >()
    );
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

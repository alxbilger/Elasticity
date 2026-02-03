#define ELASTICITY_COMPONENTS_ELEMENT_MASS_CPP

#include <Elasticity/component/ElementMass.inl>

#include <Elasticity/finiteelement/FiniteElement[all].h>
#include <sofa/core/ObjectFactory.h>

namespace elasticity
{

void registerElementMass(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData("Mass on edges")
        .add< ElementMass<sofa::defaulttype::Vec1Types, sofa::geometry::Edge> >()
        .add< ElementMass<sofa::defaulttype::Vec2Types, sofa::geometry::Edge> >()
        .add< ElementMass<sofa::defaulttype::Vec3Types, sofa::geometry::Edge> >(true));

    factory->registerObjects(sofa::core::ObjectRegistrationData("Mass on triangles")
        .add< ElementMass<sofa::defaulttype::Vec2Types, sofa::geometry::Triangle> >()
        .add< ElementMass<sofa::defaulttype::Vec3Types, sofa::geometry::Triangle> >(true));

    factory->registerObjects(sofa::core::ObjectRegistrationData("Mass on quadrangles")
        .add< ElementMass<sofa::defaulttype::Vec2Types, sofa::geometry::Quad> >()
        .add< ElementMass<sofa::defaulttype::Vec3Types, sofa::geometry::Quad> >(true));

    factory->registerObjects(sofa::core::ObjectRegistrationData("Mass on tetrahedra")
        .add< ElementMass<sofa::defaulttype::Vec3Types, sofa::geometry::Tetrahedron> >(true));

    factory->registerObjects(sofa::core::ObjectRegistrationData("Mass on hexahedra")
        .add< ElementMass<sofa::defaulttype::Vec3Types, sofa::geometry::Hexahedron> >(true));
}

template class ELASTICITY_API ElementMass<sofa::defaulttype::Vec1Types, sofa::geometry::Edge>;
template class ELASTICITY_API ElementMass<sofa::defaulttype::Vec2Types, sofa::geometry::Edge>;
template class ELASTICITY_API ElementMass<sofa::defaulttype::Vec3Types, sofa::geometry::Edge>;
template class ELASTICITY_API ElementMass<sofa::defaulttype::Vec2Types, sofa::geometry::Triangle>;
template class ELASTICITY_API ElementMass<sofa::defaulttype::Vec3Types, sofa::geometry::Triangle>;
template class ELASTICITY_API ElementMass<sofa::defaulttype::Vec2Types, sofa::geometry::Quad>;
template class ELASTICITY_API ElementMass<sofa::defaulttype::Vec3Types, sofa::geometry::Quad>;
template class ELASTICITY_API ElementMass<sofa::defaulttype::Vec3Types, sofa::geometry::Tetrahedron>;
template class ELASTICITY_API ElementMass<sofa::defaulttype::Vec3Types, sofa::geometry::Hexahedron>;

}

#include <Elasticity/LinearSmallStrainFEMForceField.inl>

#include <Elasticity/FiniteElement[Edge].h>
#include <Elasticity/FiniteElement[Triangle].h>
#include <Elasticity/FiniteElement[Quad].h>
#include <Elasticity/FiniteElement[Tetrahedron].h>
#include <Elasticity/FiniteElement[Hexahedron].h>

namespace elasticity
{

void registerLinearSmallStrainFEMForceField(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData("Hooke's law on linear beams assuming small strain")
        .add< LinearSmallStrainFEMForceField<sofa::defaulttype::Vec1Types, sofa::geometry::Edge> >()
        .add< LinearSmallStrainFEMForceField<sofa::defaulttype::Vec2Types, sofa::geometry::Edge> >()
        .add< LinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Edge> >(true));

    factory->registerObjects(sofa::core::ObjectRegistrationData("Hooke's law on linear triangles assuming small strain")
        .add< LinearSmallStrainFEMForceField<sofa::defaulttype::Vec2Types, sofa::geometry::Triangle> >()
        .add< LinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Triangle> >(true));

    factory->registerObjects(sofa::core::ObjectRegistrationData("Hooke's law on linear quads assuming small strain")
        .add< LinearSmallStrainFEMForceField<sofa::defaulttype::Vec2Types, sofa::geometry::Quad> >()
        .add< LinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Quad> >(true));

    factory->registerObjects(sofa::core::ObjectRegistrationData("Hooke's law on linear tetrahedra assuming small strain")
        .add< LinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Tetrahedron> >(true));

    factory->registerObjects(sofa::core::ObjectRegistrationData("Hooke's law on linear hexahedra assuming small strain")
        .add< LinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Hexahedron> >(true));
}

template class ELASTICITY_API LinearSmallStrainFEMForceField<sofa::defaulttype::Vec1Types, sofa::geometry::Edge>;
template class ELASTICITY_API LinearSmallStrainFEMForceField<sofa::defaulttype::Vec2Types, sofa::geometry::Edge>;
template class ELASTICITY_API LinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Edge>;
template class ELASTICITY_API LinearSmallStrainFEMForceField<sofa::defaulttype::Vec2Types, sofa::geometry::Triangle>;
template class ELASTICITY_API LinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Triangle>;
template class ELASTICITY_API LinearSmallStrainFEMForceField<sofa::defaulttype::Vec2Types, sofa::geometry::Quad>;
template class ELASTICITY_API LinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Quad>;
template class ELASTICITY_API LinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Tetrahedron>;
template class ELASTICITY_API LinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Hexahedron>;

}

#include <Elasticity/impl/CorotationalFEM.inl>

namespace elasticity
{

template class ELASTICITY_API CorotationalFEM<sofa::defaulttype::Vec2Types, sofa::geometry::Triangle>;
template class ELASTICITY_API CorotationalFEM<sofa::defaulttype::Vec2Types, sofa::geometry::Quad>;

template class ELASTICITY_API CorotationalFEM<sofa::defaulttype::Vec3Types, sofa::geometry::Triangle>;
template class ELASTICITY_API CorotationalFEM<sofa::defaulttype::Vec3Types, sofa::geometry::Quad>;
template class ELASTICITY_API CorotationalFEM<sofa::defaulttype::Vec3Types, sofa::geometry::Tetrahedron>;
template class ELASTICITY_API CorotationalFEM<sofa::defaulttype::Vec3Types, sofa::geometry::Hexahedron>;
}

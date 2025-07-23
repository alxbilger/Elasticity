#include <Elasticity/CorotationalFEM.inl>

namespace elasticity
{
template class ELASTICITY_API CorotationalFEM<sofa::defaulttype::Vec3Types, sofa::geometry::Tetrahedron>;
template class ELASTICITY_API CorotationalFEM<sofa::defaulttype::Vec3Types, sofa::geometry::Hexahedron>;
}

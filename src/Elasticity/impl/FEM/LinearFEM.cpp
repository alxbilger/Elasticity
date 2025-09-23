#define ELASTICITY_LINEARFEM_CPP

#include <Elasticity/config.h>
#include <Elasticity/finiteelement/FiniteElement[all].h>
#include <Elasticity/impl/FEM/LinearFEM.inl>
#include <sofa/defaulttype/VecTypes.h>

namespace elasticity
{

template class ELASTICITY_API LinearFEM<sofa::defaulttype::Vec1Types, sofa::geometry::Edge>;
template class ELASTICITY_API LinearFEM<sofa::defaulttype::Vec2Types, sofa::geometry::Edge>;
template class ELASTICITY_API LinearFEM<sofa::defaulttype::Vec3Types, sofa::geometry::Edge>;
template class ELASTICITY_API LinearFEM<sofa::defaulttype::Vec2Types, sofa::geometry::Triangle>;
template class ELASTICITY_API LinearFEM<sofa::defaulttype::Vec3Types, sofa::geometry::Triangle>;
template class ELASTICITY_API LinearFEM<sofa::defaulttype::Vec2Types, sofa::geometry::Quad>;
template class ELASTICITY_API LinearFEM<sofa::defaulttype::Vec3Types, sofa::geometry::Quad>;
template class ELASTICITY_API LinearFEM<sofa::defaulttype::Vec3Types, sofa::geometry::Tetrahedron>;
template class ELASTICITY_API LinearFEM<sofa::defaulttype::Vec3Types, sofa::geometry::Hexahedron>;

}

#define ELASTICITY_LINEARFEM_CPP

#include <Elasticity/LinearFEM.inl>
#include <Elasticity/config.h>
#include <sofa/defaulttype/VecTypes.h>

#include <Elasticity/FiniteElement[Edge].h>
#include <Elasticity/FiniteElement[Hexahedron].h>
#include <Elasticity/FiniteElement[Quad].h>
#include <Elasticity/FiniteElement[Tetrahedron].h>
#include <Elasticity/FiniteElement[Triangle].h>

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

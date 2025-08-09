#pragma once

#include <Elasticity/config.h>
#include <Elasticity/impl/LinearFEM.h>
#include <Elasticity/finiteelement/FiniteElement[all].h>

namespace elasticity
{
#if !defined(ELASTICITY_LINEARFEM_CPP)
extern template class ELASTICITY_API LinearFEM<sofa::defaulttype::Vec1Types, sofa::geometry::Edge>;
extern template class ELASTICITY_API LinearFEM<sofa::defaulttype::Vec2Types, sofa::geometry::Edge>;
extern template class ELASTICITY_API LinearFEM<sofa::defaulttype::Vec3Types, sofa::geometry::Edge>;
extern template class ELASTICITY_API LinearFEM<sofa::defaulttype::Vec2Types, sofa::geometry::Triangle>;
extern template class ELASTICITY_API LinearFEM<sofa::defaulttype::Vec3Types, sofa::geometry::Triangle>;
extern template class ELASTICITY_API LinearFEM<sofa::defaulttype::Vec2Types, sofa::geometry::Quad>;
extern template class ELASTICITY_API LinearFEM<sofa::defaulttype::Vec3Types, sofa::geometry::Quad>;
extern template class ELASTICITY_API LinearFEM<sofa::defaulttype::Vec3Types, sofa::geometry::Tetrahedron>;
extern template class ELASTICITY_API LinearFEM<sofa::defaulttype::Vec3Types, sofa::geometry::Hexahedron>;
#endif
}

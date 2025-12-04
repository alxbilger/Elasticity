#define ELASTICITY_COMPONENT_BASE_ELEMENT_LINEAR_FEM_FORCEFIELD_CPP
#include <Elasticity/component/BaseElementLinearFEMForceField.inl>
#include <Elasticity/finiteelement/FiniteElement[all].h>

namespace elasticity
{
template class ELASTICITY_API BaseElementLinearFEMForceField<sofa::defaulttype::Vec1Types, sofa::geometry::Edge>;
template class ELASTICITY_API BaseElementLinearFEMForceField<sofa::defaulttype::Vec2Types, sofa::geometry::Edge>;
template class ELASTICITY_API BaseElementLinearFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Edge>;
template class ELASTICITY_API BaseElementLinearFEMForceField<sofa::defaulttype::Vec2Types, sofa::geometry::Triangle>;
template class ELASTICITY_API BaseElementLinearFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Triangle>;
template class ELASTICITY_API BaseElementLinearFEMForceField<sofa::defaulttype::Vec2Types, sofa::geometry::Quad>;
template class ELASTICITY_API BaseElementLinearFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Quad>;
template class ELASTICITY_API BaseElementLinearFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Tetrahedron>;
template class ELASTICITY_API BaseElementLinearFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Hexahedron>;
}

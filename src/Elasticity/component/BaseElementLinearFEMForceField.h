#pragma once

#include <Elasticity/component/LinearMechanicalParametersComponent.h>
#include <Elasticity/config.h>
#include <Elasticity/finiteelement/FiniteElement[all].h>
#include <Elasticity/impl/ComputeStrategy.h>
#include <Elasticity/impl/trait.h>
#include <sofa/component/solidmechanics/fem/elastic/BaseLinearElasticityFEMForceField.h>

namespace elasticity
{

template <class DataTypes, class ElementType>
class BaseElementLinearFEMForceField : public sofa::component::solidmechanics::fem::elastic::BaseLinearElasticityFEMForceField<DataTypes>
{
public:
    SOFA_CLASS(
        SOFA_TEMPLATE2(BaseElementLinearFEMForceField, DataTypes, ElementType),
        sofa::component::solidmechanics::fem::elastic::BaseLinearElasticityFEMForceField<DataTypes>);

    void init() override;

private:
    using trait = elasticity::trait<DataTypes, ElementType>;
    using ElementStiffness = typename trait::ElementStiffness;
    using ElasticityTensor = typename trait::ElasticityTensor;
    using StrainDisplacement = typename trait::StrainDisplacement;

protected:

    BaseElementLinearFEMForceField();

    /**
     * With linear small strain, the element stiffness matrix is constant, so it can be precomputed.
     */
    void precomputeElementStiffness();

public:

    /**
     * List of precomputed element stiffness matrices
     */
    sofa::Data<sofa::type::vector<ElementStiffness> > d_elementStiffness;
};

#if !defined(ELASTICITY_COMPONENT_BASE_ELEMENT_LINEAR_FEM_FORCEFIELD_CPP)
extern template class ELASTICITY_API BaseElementLinearFEMForceField<sofa::defaulttype::Vec1Types, sofa::geometry::Edge>;
extern template class ELASTICITY_API BaseElementLinearFEMForceField<sofa::defaulttype::Vec2Types, sofa::geometry::Edge>;
extern template class ELASTICITY_API BaseElementLinearFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Edge>;
extern template class ELASTICITY_API BaseElementLinearFEMForceField<sofa::defaulttype::Vec2Types, sofa::geometry::Triangle>;
extern template class ELASTICITY_API BaseElementLinearFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Triangle>;
extern template class ELASTICITY_API BaseElementLinearFEMForceField<sofa::defaulttype::Vec2Types, sofa::geometry::Quad>;
extern template class ELASTICITY_API BaseElementLinearFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Quad>;
extern template class ELASTICITY_API BaseElementLinearFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Tetrahedron>;
extern template class ELASTICITY_API BaseElementLinearFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Hexahedron>;
#endif

}  // namespace elasticity

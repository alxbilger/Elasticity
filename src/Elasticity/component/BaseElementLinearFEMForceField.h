#pragma once

#include <Elasticity/component/LinearMechanicalParametersComponent.h>
#include <Elasticity/component/TopologyAccessor.h>
#include <Elasticity/config.h>
#include <Elasticity/finiteelement/FiniteElement[all].h>
#include <Elasticity/impl/ComputeStrategy.h>
#include <Elasticity/impl/trait.h>
#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/behavior/SingleStateAccessor.h>

namespace elasticity
{

template <class DataTypes, class ElementType>
class BaseElementLinearFEMForceField :
    public virtual TopologyAccessor,
    public LinearMechanicalParametersComponent<DataTypes>,
    public virtual sofa::core::behavior::SingleStateAccessor<DataTypes>
{
public:
    SOFA_CLASS3(
    SOFA_TEMPLATE2(BaseElementLinearFEMForceField, DataTypes, ElementType),
        TopologyAccessor,
        LinearMechanicalParametersComponent<DataTypes>,
        sofa::core::behavior::SingleStateAccessor<DataTypes>);

    void init() override;




private:
    using trait = elasticity::trait<DataTypes, ElementType>;
    using ElementStiffness = typename trait::ElementStiffness;
    using ElasticityTensor = typename trait::ElasticityTensor;
    using StrainDisplacement = typename trait::StrainDisplacement;
    using TopologyAccessor::l_topology;

public:
    const sofa::type::vector<ElementStiffness>& getElementStiffness() const { return m_elementStiffness; }

protected:

    BaseElementLinearFEMForceField();

    /**
     * With linear small strain, the element stiffness matrix is constant, so it can be precomputed.
     */
    void precomputeElementStiffness();

    /**
     * List of precomputed element stiffness matrices
     */
    sofa::type::vector<ElementStiffness> m_elementStiffness;

    ElasticityTensor m_elasticityTensor;

    sofa::type::vector<std::array<StrainDisplacement, trait::NumberOfNodesInElement>> m_strainDisplacement;
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

#pragma once

#include <Elasticity/component/BaseElementLinearFEMForceField.h>
#include <Elasticity/impl/ComputeStrategy.h>
#include <sofa/core/behavior/ForceField.h>
#include <sofa/helper/OptionsGroup.h>
#include <sofa/helper/decompose.h>

#include <Elasticity/component/FEMForceField.h>

#if !defined(ELASTICITY_COMPONENT_ELEMENT_COROTATIONAL_FEM_FORCE_FIELD_CPP)
#include <Elasticity/finiteelement/FiniteElement[all].h>
#endif

namespace elasticity
{

template <class DataTypes, class ElementType>
class ElementCorotationalFEMForceField :
    public BaseElementLinearFEMForceField<DataTypes, ElementType>,
    public FEMForceField<DataTypes, ElementType>
{
public:
    SOFA_CLASS2(
        SOFA_TEMPLATE2(ElementCorotationalFEMForceField, DataTypes, ElementType),
        SOFA_TEMPLATE2(BaseElementLinearFEMForceField, DataTypes, ElementType),
        SOFA_TEMPLATE2(FEMForceField, DataTypes, ElementType));

    /**
     * The purpose of this function is to register the name of this class according to the provided
     * pattern.
     *
     * Example: ElementCorotationalFEMForceField<Vec3Types, sofa::geometry::Edge> will produce
     * the class name "EdgeCorotationalFEMForceField".
     */
    static const std::string GetCustomClassName()
    {
        return std::string(sofa::geometry::elementTypeToString(ElementType::Element_type)) +
               "CorotationalFEMForceField";
    }

    static const std::string GetCustomTemplateName() { return DataTypes::Name(); }

private:
    using trait = elasticity::trait<DataTypes, ElementType>;
    using ElementForce = trait::ElementForce;
    using RotationMatrix = sofa::type::Mat<trait::spatial_dimensions, trait::spatial_dimensions, sofa::Real_t<DataTypes>>;


public:

    void init() override;

    void buildStiffnessMatrix(sofa::core::behavior::StiffnessMatrix* matrix) override;

    SReal getPotentialEnergy(const sofa::core::MechanicalParams*, const sofa::DataVecCoord_t<DataTypes>& x) const override;

    const sofa::type::vector<RotationMatrix>& getElementRotations() const { return m_rotations; }

protected:

    template<class ExecutionPolicy>
    void computeElementForce(
        sofa::type::vector<ElementForce>& elementForces,
        const sofa::VecCoord_t<DataTypes>& nodePositions);

    template<class ExecutionPolicy>
    void computeElementForceDeriv(
        sofa::type::vector<ElementForce>& elementForcesDeriv,
        const sofa::VecCoord_t<DataTypes>& nodeDx,
        sofa::Real_t<DataTypes> kFactor);

    sofa::type::vector<RotationMatrix> m_rotations;

    void computeElementRotation(
        const std::array<sofa::Coord_t<DataTypes>, trait::NumberOfNodesInElement>& nodesPosition,
        const std::array<sofa::Coord_t<DataTypes>, trait::NumberOfNodesInElement>& nodesRestPosition,
        RotationMatrix& rotationMatrix);

    sofa::Coord_t<DataTypes> translation(const std::array<sofa::Coord_t<DataTypes>, trait::NumberOfNodesInElement>& nodes) const;
    static sofa::Coord_t<DataTypes> computeCentroid(const std::array<sofa::Coord_t<DataTypes>, trait::NumberOfNodesInElement>& nodes);


};



#if !defined(ELASTICITY_COMPONENT_ELEMENT_COROTATIONAL_FEM_FORCE_FIELD_CPP)
// extern template class ELASTICITY_API ElementCorotationalFEMForceField<sofa::defaulttype::Vec1Types, sofa::geometry::Edge>;
// extern template class ELASTICITY_API ElementCorotationalFEMForceField<sofa::defaulttype::Vec2Types, sofa::geometry::Edge>;
// extern template class ELASTICITY_API ElementCorotationalFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Edge>;
extern template class ELASTICITY_API ElementCorotationalFEMForceField<sofa::defaulttype::Vec2Types, sofa::geometry::Triangle>;
extern template class ELASTICITY_API ElementCorotationalFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Triangle>;
extern template class ELASTICITY_API ElementCorotationalFEMForceField<sofa::defaulttype::Vec2Types, sofa::geometry::Quad>;
extern template class ELASTICITY_API ElementCorotationalFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Quad>;
extern template class ELASTICITY_API ElementCorotationalFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Tetrahedron>;
extern template class ELASTICITY_API ElementCorotationalFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Hexahedron>;
#endif

}  // namespace elasticity

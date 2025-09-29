#pragma once

#include <Elasticity/component/ElementLinearSmallStrainFEMForceField.h>
#include <sofa/helper/decompose.h>

#if !defined(ELASTICITY_COMPONENT_ELEMENT_COROTATIONAL_FEM_FORCE_FIELD_CPP)
#include <Elasticity/finiteelement/FiniteElement[all].h>
#endif

namespace elasticity
{

template <class DataTypes, class ElementType>
class ElementCorotationalFEMForceField : public ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>
{
public:
    SOFA_CLASS(
        SOFA_TEMPLATE2(ElementCorotationalFEMForceField, DataTypes, ElementType),
        SOFA_TEMPLATE2(ElementLinearSmallStrainFEMForceField, DataTypes, ElementType));

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
    using DataVecCoord = sofa::DataVecDeriv_t<DataTypes>;
    using DataVecDeriv = sofa::DataVecDeriv_t<DataTypes>;
    using VecCoord = sofa::VecCoord_t<DataTypes>;
    using VecDeriv = sofa::VecDeriv_t<DataTypes>;
    using Coord = sofa::Coord_t<DataTypes>;
    using Deriv = sofa::Deriv_t<DataTypes>;
    using Real = sofa::Real_t<DataTypes>;

    using FiniteElement = elasticity::FiniteElement<ElementType, DataTypes>;
    using ReferenceCoord = typename FiniteElement::ReferenceCoord;

    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;
    static constexpr sofa::Size NumberOfNodesInElement = ElementType::NumberOfNodes;
    static constexpr sofa::Size NumberOfDofsInElement = NumberOfNodesInElement * spatial_dimensions;
    static constexpr sofa::Size ElementDimension = FiniteElement::ElementDimension;

    /// the concatenation of the displacement of the element nodes in a single vector
    using ElementDisplacement = sofa::type::Vec<NumberOfDofsInElement, Real>;

    /// the type of the element stiffness matrix
    using ElementStiffness = ElementStiffness<DataTypes, ElementType>;

    using TopologyAccessor::l_topology;

public:

    void addForce(const sofa::core::MechanicalParams* mparams, DataVecDeriv& f,
              const DataVecCoord& x, const DataVecDeriv& v) override;

    void addDForce(const sofa::core::MechanicalParams* mparams, DataVecDeriv& df,
                   const DataVecDeriv& dx) override;

    void buildStiffnessMatrix(sofa::core::behavior::StiffnessMatrix* matrix) override;

    SReal getPotentialEnergy(const sofa::core::MechanicalParams*, const DataVecCoord& x) const override;

protected:

    using RotationMatrix = sofa::type::Mat<spatial_dimensions, spatial_dimensions, Real>;
    sofa::type::vector<RotationMatrix> m_rotations;

    void computeElementRotation(
        const std::array<Coord, NumberOfNodesInElement>& nodesPosition,
        const std::array<Coord, NumberOfNodesInElement>& nodesRestPosition,
        RotationMatrix& rotationMatrix);

    Coord translation(const std::array<Coord, NumberOfNodesInElement>& nodes) const;
    static Coord computeCentroid(const std::array<Coord, NumberOfNodesInElement>& nodes);
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

#pragma once

#include <Elasticity/config.h>
#include <Elasticity/impl/ElementStiffnessMatrix.h>
#include <sofa/core/behavior/ForceField.h>

#include <Elasticity/component/LinearMechanicalParametersComponent.h>
#include <Elasticity/component/TopologyAccessor.h>

#if !defined(ELASTICITY_COMPONENT_ELEMENT_LINEAR_SMALL_STRAIN_FEM_FORCE_FIELD_CPP)
#include <Elasticity/finiteelement/FiniteElement[all].h>
#endif

namespace elasticity
{

template <class DataTypes, class ElementType>
class ElementLinearSmallStrainFEMForceField :
    public TopologyAccessor,
    public LinearMechanicalParametersComponent<DataTypes>,
    public sofa::core::behavior::ForceField<DataTypes>
{
public:
    SOFA_CLASS3(
        SOFA_TEMPLATE2(ElementLinearSmallStrainFEMForceField, DataTypes, ElementType),
            TopologyAccessor,
            LinearMechanicalParametersComponent<DataTypes>,
            sofa::core::behavior::ForceField<DataTypes>);

    /**
     * The purpose of this function is to register the name of this class according to the provided
     * pattern.
     *
     * Example: ElementLinearSmallStrainFEMForceField<Vec3Types, sofa::geometry::Edge> will produce
     * the class name "EdgeLinearSmallStrainFEMForceField".
     */
    static const std::string GetCustomClassName()
    {
        return std::string(sofa::geometry::elementTypeToString(ElementType::Element_type)) +
               "LinearSmallStrainFEMForceField";
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

    /// type of 2nd-order tensor for the elasticity tensor for isotropic materials
    using ElasticityTensor = elasticity::ElasticityTensor<DataTypes>;

    /// the type of B in e = B d, if e is the strain, and d is the displacement
    using StrainDisplacement = elasticity::StrainDisplacement<DataTypes, ElementType>;

    /// the concatenation of the displacement of the element nodes in a single vector
    using ElementDisplacement = sofa::type::Vec<NumberOfDofsInElement, Real>;

    /// the type of the element stiffness matrix
    using ElementStiffness = elasticity::ElementStiffness<DataTypes, ElementType>;

public:
    void init() override;

    void addForce(const sofa::core::MechanicalParams* mparams, DataVecDeriv& f,
              const DataVecCoord& x, const DataVecDeriv& v) override;

    void addDForce(const sofa::core::MechanicalParams* mparams, DataVecDeriv& df,
                   const DataVecDeriv& dx) override;

    void buildStiffnessMatrix(sofa::core::behavior::StiffnessMatrix* matrix) override;

    SReal getPotentialEnergy(const sofa::core::MechanicalParams*, const DataVecCoord& x) const override;

    using sofa::core::behavior::ForceField<DataTypes>::addKToMatrix;
    // almost deprecated, but here for compatibility with unit tests
    void addKToMatrix(sofa::linearalgebra::BaseMatrix* matrix, SReal kFact, unsigned& offset) override;

   protected:

    /**
     * With linear small strain, the element stiffness matrix is constant, so it can be precomputed.
     */
    void precomputeElementStiffness();

    /**
     * List of precomputed element stiffness matrices
     */
    sofa::type::vector<ElementStiffness> m_elementStiffness;

    ElasticityTensor m_elasticityTensor;

    sofa::type::vector<std::array<StrainDisplacement, NumberOfNodesInElement>> m_strainDisplacement;

    static ElementDisplacement computeElementDisplacement(
        const std::array<Coord, NumberOfNodesInElement>& elementNodesCoordinates,
        const std::array<Coord, NumberOfNodesInElement>& restElementNodesCoordinates);
};

#if !defined(ELASTICITY_COMPONENT_ELEMENT_LINEAR_SMALL_STRAIN_FEM_FORCE_FIELD_CPP)
extern template class ELASTICITY_API ElementLinearSmallStrainFEMForceField<sofa::defaulttype::Vec1Types, sofa::geometry::Edge>;
extern template class ELASTICITY_API ElementLinearSmallStrainFEMForceField<sofa::defaulttype::Vec2Types, sofa::geometry::Edge>;
extern template class ELASTICITY_API ElementLinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Edge>;
extern template class ELASTICITY_API ElementLinearSmallStrainFEMForceField<sofa::defaulttype::Vec2Types, sofa::geometry::Triangle>;
extern template class ELASTICITY_API ElementLinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Triangle>;
extern template class ELASTICITY_API ElementLinearSmallStrainFEMForceField<sofa::defaulttype::Vec2Types, sofa::geometry::Quad>;
extern template class ELASTICITY_API ElementLinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Quad>;
extern template class ELASTICITY_API ElementLinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Tetrahedron>;
extern template class ELASTICITY_API ElementLinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Hexahedron>;
#endif

}  // namespace elasticity

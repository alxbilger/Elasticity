#pragma once

#include <Elasticity/component/HyperelasticMaterial.h>
#include <Elasticity/config.h>
#include <Elasticity/finiteelement/FiniteElement.h>
#include <Elasticity/impl/ElementStiffnessMatrix.h>
#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/behavior/TopologyAccessor.h>

#if !defined(ELASTICITY_COMPONENT_ELEMENT_HYPERLASTICITY_FEM_FORCE_FIELD_CPP)
#include <Elasticity/finiteelement/FiniteElement[all].h>
#endif

namespace elasticity
{

template <class TDataTypes, class TElementType>
class ElementHyperelasticityFEMForceField :
    public sofa::core::behavior::TopologyAccessor,
    public sofa::core::behavior::ForceField<TDataTypes>
{
public:
    SOFA_CLASS2(
        SOFA_TEMPLATE2(ElementHyperelasticityFEMForceField, TDataTypes, TElementType),
            sofa::core::behavior::TopologyAccessor,
            sofa::core::behavior::ForceField<TDataTypes>);

    using DataTypes = TDataTypes;

    /**
     * The purpose of this function is to register the name of this class according to the provided
     * pattern.
     *
     * Example: ElementHyperelasticityFEMForceField<Vec3Types, sofa::geometry::Edge> will produce
     * the class name "EdgeHyperelasticityFEMForceField".
     */
    static const std::string GetCustomClassName()
    {
        return std::string(sofa::geometry::elementTypeToString(TElementType::Element_type)) +
               "HyperelasticityFEMForceField";
    }

    static const std::string GetCustomTemplateName() { return TDataTypes::Name(); }

private:
    using DataVecCoord = sofa::DataVecDeriv_t<TDataTypes>;
    using DataVecDeriv = sofa::DataVecDeriv_t<TDataTypes>;
    using VecCoord = sofa::VecCoord_t<TDataTypes>;
    using VecDeriv = sofa::VecDeriv_t<TDataTypes>;
    using Coord = sofa::Coord_t<TDataTypes>;
    using Deriv = sofa::Deriv_t<TDataTypes>;
    using Real = sofa::Real_t<TDataTypes>;

    using FiniteElement = elasticity::FiniteElement<TElementType, TDataTypes>;

    static constexpr sofa::Size spatial_dimensions = TDataTypes::spatial_dimensions;
    static constexpr sofa::Size NumberOfNodesInElement = TElementType::NumberOfNodes;
    static constexpr sofa::Size NumberOfDofsInElement = NumberOfNodesInElement * spatial_dimensions;
    static constexpr sofa::Size TopologicalDimension = FiniteElement::TopologicalDimension;
    static constexpr sofa::Size NumberOfQuadraturePoints = FiniteElement::quadraturePoints().size();

    using DeformationGradient = sofa::type::Mat<spatial_dimensions, spatial_dimensions, Real>;

    /// the type of the element stiffness matrix
    using ElementStiffness = elasticity::ElementStiffness<TDataTypes, TElementType>;

public:
    void init() override;

    void addForce(const sofa::core::MechanicalParams* mparams, DataVecDeriv& f,
          const DataVecCoord& x, const DataVecDeriv& v) override;

    void addDForce(const sofa::core::MechanicalParams* mparams, DataVecDeriv& df,
                   const DataVecDeriv& dx) override;

    void buildStiffnessMatrix(sofa::core::behavior::StiffnessMatrix* matrix) override;

    using sofa::core::behavior::ForceField<DataTypes>::getPotentialEnergy;
    SReal getPotentialEnergy(const sofa::core::MechanicalParams*, const DataVecCoord& x) const override;

    using sofa::core::behavior::ForceField<DataTypes>::addKToMatrix;
    // almost deprecated, but here for compatibility with unit tests
    void addKToMatrix(sofa::linearalgebra::BaseMatrix* matrix, SReal kFact, unsigned& offset) override;

    sofa::SingleLink<MyType, HyperelasticMaterial<TDataTypes>,
        sofa::BaseLink::FLAG_STOREPATH | sofa::BaseLink::FLAG_STRONGLINK> l_material;

protected:
    void validateMaterial();

    bool m_isHessianValid;

    void computeHessian(const VecCoord& coordinates);

    /**
     * List of precomputed element stiffness matrices
     */
    sofa::type::vector<ElementStiffness> m_elementStiffness;

    const VecCoord* m_coordinates{ nullptr };

    DeformationGradient computeDeformationGradient(
        const sofa::type::Mat<spatial_dimensions, TopologicalDimension, Real>& J_q,
        const sofa::type::Mat<TopologicalDimension, spatial_dimensions, Real>& J_Q_inv);
    DeformationGradient computeDeformationGradient2(
        const std::array<Coord, NumberOfNodesInElement>& elementNodesCoordinates,
        const sofa::type::Mat<NumberOfNodesInElement, spatial_dimensions, Real>& dN_dQ);

    struct PrecomputedData
    {
        // jacobian of the mapping from the reference space to the rest physical space, evaluated at the
        // quadrature point
        sofa::type::Mat<spatial_dimensions, TopologicalDimension, Real> jacobian { sofa::type::NOINIT };

        // inverse of the jacobian of the mapping from the reference space to the rest physical space,
        // evaluated at the quadrature point
        sofa::type::Mat<TopologicalDimension, spatial_dimensions, Real> jacobianInv { sofa::type::NOINIT };

        Real detJacobian {};

        // gradient of the shape functions in the physical element evaluated at the quadrature point
        sofa::type::Mat<NumberOfNodesInElement, spatial_dimensions, Real> dN_dQ { sofa::type::NOINIT };
    };

    /**
     * Data can be precomputed from the rest configuration and used later in computations for the
     * current configuration. A piece of precomputed data is stored for each quadrature point in
     * each element.
     */
    sofa::type::vector<std::array<PrecomputedData, NumberOfQuadraturePoints>> m_precomputedData;

    void precomputeData();
};

#if !defined(ELASTICITY_COMPONENT_ELEMENT_HYPERLASTICITY_FEM_FORCE_FIELD_CPP)
template class ELASTICITY_API ElementHyperelasticityFEMForceField<sofa::defaulttype::Vec1Types, sofa::geometry::Edge>;
template class ELASTICITY_API ElementHyperelasticityFEMForceField<sofa::defaulttype::Vec2Types, sofa::geometry::Edge>;
template class ELASTICITY_API ElementHyperelasticityFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Edge>;
template class ELASTICITY_API ElementHyperelasticityFEMForceField<sofa::defaulttype::Vec2Types, sofa::geometry::Triangle>;
template class ELASTICITY_API ElementHyperelasticityFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Triangle>;
template class ELASTICITY_API ElementHyperelasticityFEMForceField<sofa::defaulttype::Vec2Types, sofa::geometry::Quad>;
template class ELASTICITY_API ElementHyperelasticityFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Quad>;
template class ELASTICITY_API ElementHyperelasticityFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Tetrahedron>;
template class ELASTICITY_API ElementHyperelasticityFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Hexahedron>;
#endif

}  // namespace elasticity

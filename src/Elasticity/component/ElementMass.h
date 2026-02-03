#pragma once

#include <Elasticity/config.h>
#include <Elasticity/impl/trait.h>
#include <sofa/core/behavior/Mass.h>
#include <sofa/core/behavior/TopologyAccessor.h>

#if !defined(ELASTICITY_COMPONENTS_ELEMENT_MASS_CPP)
#include <Elasticity/finiteelement/FiniteElement[all].h>
#endif

namespace elasticity
{

template <class DataTypes, class ElementType>
class ElementMass :
    public sofa::core::behavior::Mass<DataTypes>,
    public sofa::core::behavior::TopologyAccessor
{
public:
    SOFA_CLASS2(SOFA_TEMPLATE2(ElementMass, DataTypes, ElementType),
           sofa::core::behavior::Mass<DataTypes>,
           sofa::core::behavior::TopologyAccessor);

    /**
     * The purpose of this function is to register the name of this class according to the provided
     * pattern.
     *
     * Example: ElementMass<Vec3Types, sofa::geometry::Edge> will produce
     * the class name "EdgeMass".
     */
    static const std::string GetCustomClassName()
    {
        return std::string(sofa::geometry::elementTypeToString(ElementType::Element_type)) +
               "Mass";
    }

    static const std::string GetCustomTemplateName() { return DataTypes::Name(); }

    bool isDiagonal() const override;

    void init() override;

    void addForce(
        const sofa::core::MechanicalParams* mparams,
        sofa::DataVecDeriv_t<DataTypes>& f,
        const sofa::DataVecCoord_t<DataTypes>& x,
        const sofa::DataVecDeriv_t<DataTypes>& v) override;

    void addMDx(const sofa::core::MechanicalParams* mparams,
        sofa::DataVecDeriv_t<DataTypes>& vres,
        const sofa::DataVecDeriv_t<DataTypes>& vdx,
        SReal factor) override;

    void addMToMatrix(sofa::linearalgebra::BaseMatrix * matrix, SReal mFact, unsigned int &offset) override;

    sofa::Data<sofa::VecReal_t<DataTypes> > d_nodalDensity;

protected:

    using trait = elasticity::trait<DataTypes, ElementType>;
    using FiniteElement = typename trait::FiniteElement;

    using Coord = sofa::Coord_t<DataTypes>;
    using Deriv = sofa::Deriv_t<DataTypes>;
    using Real = sofa::Real_t<DataTypes>;

    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;
    static constexpr sofa::Size NumberOfNodesInElement = trait::NumberOfNodesInElement;

    using ScalarMassMatrix = sofa::type::Mat<NumberOfNodesInElement, NumberOfNodesInElement, Real>;
    using ElementQuadratureMass = std::array<ScalarMassMatrix, trait::NbQuadraturePoints>;

    ElementMass();

    static constexpr Real defaultNodalDensity { 1. };
    void resizeNodalDensity(const std::size_t size);

    Deriv getGravity() const;

    void updateMassCacheIfNeeded();

    sofa::type::vector<ElementQuadratureMass> m_elementQuadratureMass; // geometry-only cache (rest positions)
};



#if !defined(ELASTICITY_COMPONENTS_ELEMENT_MASS_CPP)
extern template class ELASTICITY_API ElementMass<sofa::defaulttype::Vec1Types, sofa::geometry::Edge>;
extern template class ELASTICITY_API ElementMass<sofa::defaulttype::Vec2Types, sofa::geometry::Edge>;
extern template class ELASTICITY_API ElementMass<sofa::defaulttype::Vec3Types, sofa::geometry::Edge>;
extern template class ELASTICITY_API ElementMass<sofa::defaulttype::Vec2Types, sofa::geometry::Triangle>;
extern template class ELASTICITY_API ElementMass<sofa::defaulttype::Vec3Types, sofa::geometry::Triangle>;
extern template class ELASTICITY_API ElementMass<sofa::defaulttype::Vec2Types, sofa::geometry::Quad>;
extern template class ELASTICITY_API ElementMass<sofa::defaulttype::Vec3Types, sofa::geometry::Quad>;
extern template class ELASTICITY_API ElementMass<sofa::defaulttype::Vec3Types, sofa::geometry::Tetrahedron>;
extern template class ELASTICITY_API ElementMass<sofa::defaulttype::Vec3Types, sofa::geometry::Hexahedron>;
#endif

}

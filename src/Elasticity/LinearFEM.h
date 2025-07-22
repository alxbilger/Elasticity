#pragma once

#include <Elasticity/config.h>
#include <Elasticity/FiniteElement.h>
#include <sofa/core/behavior/BaseLocalForceFieldMatrix.h>

namespace elasticity
{

template <class DataTypes>
class BaseLinearFEM
{
    using VecCoord = sofa::VecCoord_t<DataTypes>;
    using VecDeriv = sofa::VecDeriv_t<DataTypes>;
    using Real = sofa::Real_t<DataTypes>;

public:
    ~BaseLinearFEM() = default;

    virtual void addForce(VecDeriv& force, const VecCoord& position, const VecCoord& restPosition) const = 0;
    virtual void addDForce(VecDeriv& df, const VecDeriv& dx, Real kFactor) const = 0;
    virtual void buildStiffnessMatrix(sofa::core::behavior::StiffnessMatrix::Derivative& dfdx) const = 0;
    virtual void precomputeElementStiffness(const VecCoord& restPosition, Real youngModulus, Real poissonRatio) = 0;
};

template <class DataTypes, class ElementType>
class LinearFEM final : public BaseLinearFEM<DataTypes>
{
    using VecCoord = sofa::VecCoord_t<DataTypes>;
    using VecDeriv = sofa::VecDeriv_t<DataTypes>;
    using Coord = sofa::Coord_t<DataTypes>;
    using Deriv = sofa::Deriv_t<DataTypes>;
    using Real = sofa::Real_t<DataTypes>;
    using TopologyElement = sofa::topology::Element<ElementType>;
    using FiniteElement = FiniteElement<ElementType, DataTypes>;

    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;
    static constexpr sofa::Size NumberOfNodesInElement = ElementType::NumberOfNodes;
    static constexpr sofa::Size NumberOfDofsInElement = NumberOfNodesInElement * spatial_dimensions;
    static constexpr sofa::Size ElementDimension = FiniteElement::ElementDimension;

    /// The number of independent elements in a symmetric 2nd-order tensor of size
    /// (spatial_dimensions x spatial_dimensions)
    static constexpr sofa::Size NumberOfIndependentElements = spatial_dimensions * (spatial_dimensions + 1) / 2;

    /// type of 2nd-order tensor for the elasticity tensor for isotropic materials
    using ElasticityTensor = sofa::type::Mat<NumberOfIndependentElements, NumberOfIndependentElements, Real>;

    /// the type of B in e = B d, if e is the strain, and d is the displacement
    using StrainDisplacement = sofa::type::Mat<NumberOfIndependentElements, NumberOfDofsInElement, Real>;

    /// the concatenation of the displacement of the 4 nodes in a single vector
    using ElementDisplacement = sofa::type::Vec<NumberOfDofsInElement, Real>;

    /// the type of the element stiffness matrix
    using ElementStiffness = sofa::type::Mat<NumberOfDofsInElement, NumberOfDofsInElement, Real>;

public:
    explicit LinearFEM(sofa::core::topology::BaseMeshTopology* topology = nullptr);

    void addForce(VecDeriv& force, const VecCoord& position, const VecCoord& restPosition) const override;
    void addDForce(VecDeriv& df, const VecDeriv& dx, Real kFactor) const override;
    void buildStiffnessMatrix(sofa::core::behavior::StiffnessMatrix::Derivative& dfdx) const override;

    void setTopology(sofa::core::topology::BaseMeshTopology* topology);

    void precomputeElementStiffness(const VecCoord& restPosition, Real youngModulus, Real poissonRatio) override;

    static ElasticityTensor computeElasticityTensor(Real youngModulus, Real poissonRatio);

    static StrainDisplacement buildStrainDisplacement(const sofa::type::Mat<NumberOfNodesInElement, spatial_dimensions, Real> gradientShapeFunctions);

protected:

    /**
     * List of precomputed element stiffness matrices
     */
    sofa::type::vector<ElementStiffness> m_elementStiffness;

    sofa::core::topology::BaseMeshTopology* m_topology { nullptr };

    /**
     * Assemble in a unique vector the displacement of all the nodes in an element
     */
    static ElementDisplacement computeElementDisplacement(
        const std::array<Coord, NumberOfNodesInElement>& elementNodesCoordinates,
        const std::array<Coord, NumberOfNodesInElement>& restElementNodesCoordinates);
};

#if !defined(ELASTICITY_LINEARFEM_CPP)
#include <Elasticity/FiniteElement[Edge].h>
#include <Elasticity/FiniteElement[Hexahedron].h>
#include <Elasticity/FiniteElement[Quad].h>
#include <Elasticity/FiniteElement[Tetrahedron].h>
#include <Elasticity/FiniteElement[Triangle].h>

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

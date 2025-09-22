#pragma once

#include <Elasticity/config.h>
#include <Elasticity/finiteelement/FiniteElement.h>
#include <Elasticity/impl/ElementStiffnessMatrix.h>
#include <Elasticity/impl/VonMisesStressContainer.h>
#include <Elasticity/impl/SymmetricTensor.h>
#include <sofa/core/behavior/BaseLocalForceFieldMatrix.h>

#include <Elasticity/impl/ElasticityTensor.h>

namespace elasticity
{

/**
 * Base class for linear FEM. The class is designed to be called in a ForceField. The derived classes
 * are specific to a type of element.
 */
template <class DataTypes>
class BaseLinearFEM
{
    using VecCoord = sofa::VecCoord_t<DataTypes>;
    using VecDeriv = sofa::VecDeriv_t<DataTypes>;
    using Real = sofa::Real_t<DataTypes>;

public:
    virtual ~BaseLinearFEM() = default;

    virtual void addForce(VecDeriv& force, const VecCoord& position, const VecCoord& restPosition) = 0;
    virtual void addDForce(VecDeriv& df, const VecDeriv& dx, Real kFactor) const = 0;
    virtual void buildStiffnessMatrix(sofa::core::behavior::StiffnessMatrix::Derivative& dfdx) const = 0;
    virtual void precomputeElementStiffness(const VecCoord& restPosition, Real youngModulus, Real poissonRatio) = 0;
    virtual void computeVonMisesStress(VonMisesStressContainer<Real>& vonMisesStressContainer, const VecCoord& position, const VecCoord& restPosition) const {}
};

template <class DataTypes, class ElementType>
class LinearFEM : public BaseLinearFEM<DataTypes>
{
protected:
    using VecCoord = sofa::VecCoord_t<DataTypes>;
    using VecDeriv = sofa::VecDeriv_t<DataTypes>;
    using Coord = sofa::Coord_t<DataTypes>;
    using Deriv = sofa::Deriv_t<DataTypes>;
    using Real = sofa::Real_t<DataTypes>;
    using TopologyElement = sofa::topology::Element<ElementType>;
    using FiniteElement = elasticity::FiniteElement<ElementType, DataTypes>;
    using ReferenceCoord = typename FiniteElement::ReferenceCoord;

    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;
    static constexpr sofa::Size NumberOfNodesInElement = ElementType::NumberOfNodes;
    static constexpr sofa::Size NumberOfDofsInElement = NumberOfNodesInElement * spatial_dimensions;
    static constexpr sofa::Size ElementDimension = FiniteElement::ElementDimension;

    /// The number of independent elements in a symmetric 2nd-order tensor of size
    /// (spatial_dimensions x spatial_dimensions)
    static constexpr sofa::Size NumberOfIndependentElements = symmetric_tensor::NumberOfIndependentElements<spatial_dimensions>;

    /// type of 2nd-order tensor for the elasticity tensor for isotropic materials
    using ElasticityTensor = ElasticityTensor<DataTypes>;

    /// the type of B in e = B d, if e is the strain, and d is the displacement
    using StrainDisplacement = StrainDisplacement<DataTypes, ElementType>;

    /// the concatenation of the displacement of the 4 nodes in a single vector
    using ElementDisplacement = sofa::type::Vec<NumberOfDofsInElement, Real>;

    /// the type of the element stiffness matrix
    using ElementStiffness = ElementStiffness<DataTypes, ElementType>;

public:
    explicit LinearFEM(sofa::core::topology::BaseMeshTopology* topology = nullptr);

    void addForce(VecDeriv& force, const VecCoord& position, const VecCoord& restPosition) override;
    void addDForce(VecDeriv& df, const VecDeriv& dx, Real kFactor) const override;
    void buildStiffnessMatrix(sofa::core::behavior::StiffnessMatrix::Derivative& dfdx) const override;

    void precomputeElementStiffness(const VecCoord& restPosition, Real youngModulus, Real poissonRatio) override;

    void computeVonMisesStress(VonMisesStressContainer<Real>& vonMisesStressContainer,
                               const VecCoord& position,
                               const VecCoord& restPosition) const override;

protected:

    /**
     * List of precomputed element stiffness matrices
     */
    sofa::type::vector<ElementStiffness> m_elementStiffness;

    sofa::type::vector<std::array<StrainDisplacement, NumberOfNodesInElement>> m_strainDisplacement;

    ElasticityTensor m_elasticityTensor;

    sofa::core::topology::BaseMeshTopology* m_topology { nullptr };

    /**
     * Assemble in a unique vector the displacement of all the nodes in an element
     */
    static ElementDisplacement computeElementDisplacement(
        const std::array<Coord, NumberOfNodesInElement>& elementNodesCoordinates,
        const std::array<Coord, NumberOfNodesInElement>& restElementNodesCoordinates);

    const sofa::type::vector<ElementStiffness>& stiffnessMatrices() const;
};

#if !defined(ELASTICITY_LINEARFEM_CPP)
#include <Elasticity/finiteelement/FiniteElement[Edge].h>
#include <Elasticity/finiteelement/FiniteElement[Hexahedron].h>
#include <Elasticity/finiteelement/FiniteElement[Quad].h>
#include <Elasticity/finiteelement/FiniteElement[Tetrahedron].h>
#include <Elasticity/finiteelement/FiniteElement[Triangle].h>

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

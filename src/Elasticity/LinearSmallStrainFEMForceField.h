#pragma once
#include <Elasticity/FiniteElement.h>
#include <Elasticity/config.h>
#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/topology/BaseMeshTopology.h>

namespace elasticity
{

template <class DataTypes, class ElementType>
class LinearSmallStrainFEMForceField : public sofa::core::behavior::ForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE2(LinearSmallStrainFEMForceField, DataTypes, ElementType),
               SOFA_TEMPLATE(sofa::core::behavior::ForceField, DataTypes));

    static const std::string GetCustomClassName()
    {
        return std::string(sofa::geometry::elementTypeToString(ElementType::Element_type)) + "LinearSmallStrainFEMForceField";
    }

    static const std::string GetCustomTemplateName()
    {
        return DataTypes::Name();
    }

private:
    using DataVecCoord = sofa::DataVecDeriv_t<DataTypes>;
    using DataVecDeriv = sofa::DataVecDeriv_t<DataTypes>;
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

    /// The number of independent elements in a symmetric 2nd-order tensor of size (spatial_dimensions x spatial_dimensions)
    static constexpr sofa::Size NumberOfIndependentElements = spatial_dimensions * (spatial_dimensions + 1) / 2;

    /// type of 2nd-order tensor for the elasticity tensor for isotropic materials
    using ElasticityTensor = sofa::type::Mat<NumberOfIndependentElements, NumberOfIndependentElements, Real>;

    /// the type of B in e = B d, if e is the strain, and d is the displacement
    using StrainDisplacement = sofa::type::Mat<NumberOfIndependentElements, NumberOfDofsInElement, Real>;

    /// the concatenation of the displacement of the 4 nodes in a single vector
    using ElementDisplacement = sofa::type::Vec<NumberOfDofsInElement, Real>;

    /// the type of the element stiffness matrix
    using ElementStiffness = sofa::type::Mat<NumberOfDofsInElement, NumberOfDofsInElement, Real>;

    LinearSmallStrainFEMForceField();

public:
    void init() override;

    void addForce(const sofa::core::MechanicalParams* mparams, DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& v) override;

    void addDForce(const sofa::core::MechanicalParams* mparams, DataVecDeriv& df, const DataVecDeriv& dx) override;

    void buildStiffnessMatrix(sofa::core::behavior::StiffnessMatrix* matrix) override;

    SReal getPotentialEnergy(const sofa::core::MechanicalParams*,
                             const DataVecCoord& x) const override;

    /// The topology will give access to the elements
    sofa::SingleLink<LinearSmallStrainFEMForceField, sofa::core::topology::BaseMeshTopology,
                     sofa::BaseLink::FLAG_STOREPATH | sofa::BaseLink::FLAG_STRONGLINK> l_topology;

    sofa::Data<Real> d_poissonRatio;
    sofa::Data<Real> d_youngModulus;

    static ElasticityTensor computeElasticityTensor(Real youngModulus, Real poissonRatio);

    static StrainDisplacement buildStrainDisplacement(const sofa::type::Mat<NumberOfNodesInElement, spatial_dimensions, Real> gradientShapeFunctions);

protected:

    ElasticityTensor computeElasticityTensor();

    /**
     * List of precomputed element stiffness matrices
     */
    sofa::type::vector<ElementStiffness> m_elementStiffness;

    /**
     * Ensure a link to a valid topology. Without a topology, this force field cannot have access
     * to the list of elements.
     */
    void validateTopology();

    /**
     * With linear small strain, the element stiffness matrix is constant, so it can be precomputed.
     */
    void precomputeElementStiffness();

    /**
     * Assemble in a unique vector the displacement of all the nodes in an element
     */
    ElementDisplacement computeElementDisplacement(
        const std::array<Coord, NumberOfNodesInElement>& elementNodesCoordinates,
        const std::array<Coord, NumberOfNodesInElement>& restElementNodesCoordinates);
};

}  // namespace elasticity

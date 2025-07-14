#pragma once
#include <Elasticity/config.h>
#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/topology/BaseMeshTopology.h>

namespace elasticity
{

template <class DataTypes>
class TET4LinearSmallStrainFEMForceField : public sofa::core::behavior::ForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(TET4LinearSmallStrainFEMForceField, DataTypes),
               SOFA_TEMPLATE(sofa::core::behavior::ForceField, DataTypes));

private:
    using DataVecCoord = sofa::DataVecDeriv_t<DataTypes>;
    using DataVecDeriv = sofa::DataVecDeriv_t<DataTypes>;
    using VecCoord = sofa::DataVecCoord_t<DataTypes>;
    using VecDeriv = sofa::DataVecDeriv_t<DataTypes>;
    using Coord = sofa::Coord_t<DataTypes>;
    using Deriv = sofa::Deriv_t<DataTypes>;
    using Real = sofa::Real_t<DataTypes>;

    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;
    static constexpr sofa::Size NumberOfNodesInElement = sofa::geometry::Tetrahedron::NumberOfNodes;
    static constexpr sofa::Size NumberOfDofsInElement = NumberOfNodesInElement * spatial_dimensions;

    using ElasticityTensor = sofa::type::Mat<6, 6, Real>;
    using StrainDisplacement = sofa::type::Mat<6, NumberOfDofsInElement, Real>;
    using ElementDisplacement = sofa::type::Vec<NumberOfDofsInElement, Real>;
    using ElementStiffness = sofa::type::Mat<NumberOfDofsInElement, NumberOfDofsInElement, Real>;

    TET4LinearSmallStrainFEMForceField();

public:
    void init() override;

    void addForce(const sofa::core::MechanicalParams* mparams, DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& v) override;

    void addDForce(const sofa::core::MechanicalParams* mparams, DataVecDeriv& df, const DataVecDeriv& dx) override;

    SReal getPotentialEnergy(const sofa::core::MechanicalParams*,
                             const DataVecCoord& x) const override;

    sofa::SingleLink<TET4LinearSmallStrainFEMForceField, sofa::core::topology::BaseMeshTopology,
                     sofa::BaseLink::FLAG_STOREPATH | sofa::BaseLink::FLAG_STRONGLINK> l_topology;

    sofa::Data<Real> d_poissonRatio;
    sofa::Data<Real> d_youngModulus;

protected:

    ElasticityTensor computeElasticityTensor();
    StrainDisplacement computeStrainDisplacement(const std::array<Coord, NumberOfNodesInElement>& tetraNodesCoordinates);

    /**
     * List of precomputed element stiffness matrices
     */
    sofa::type::vector<ElementStiffness> m_elementStiffness;

    /**
     * Ensure a link to a valid topology. Without a topology, this force field cannot have access
     * to the list of tetrahedra.
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
        const std::array<Coord, NumberOfNodesInElement>& tetraNodesCoordinates,
        const std::array<Coord, NumberOfNodesInElement>& restTetraNodesCoordinates);
};

}  // namespace elasticity

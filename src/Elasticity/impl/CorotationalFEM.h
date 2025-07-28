#pragma once

#include <Elasticity/impl/LinearFEM.h>

namespace elasticity
{

template <class DataTypes, class ElementType>
class CorotationalFEM : public LinearFEM<DataTypes, ElementType>
{
    using Real = sofa::Real_t<DataTypes>;
    using VecCoord = sofa::VecCoord_t<DataTypes>;
    using VecDeriv = sofa::VecDeriv_t<DataTypes>;
    using Coord = sofa::Coord_t<DataTypes>;
    using Deriv = sofa::Deriv_t<DataTypes>;
    using TopologyElement = typename LinearFEM<DataTypes, ElementType>::TopologyElement;
    using FiniteElement = typename LinearFEM<DataTypes, ElementType>::FiniteElement;

    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;
    static constexpr sofa::Size NumberOfNodesInElement = ElementType::NumberOfNodes;
    static constexpr sofa::Size NumberOfDofsInElement = NumberOfNodesInElement * spatial_dimensions;
    static constexpr sofa::Size ElementDimension = FiniteElement::ElementDimension;
    static constexpr sofa::Size NumberOfIndependentElements = LinearFEM<DataTypes, ElementType>::NumberOfIndependentElements;

    using ElementStiffness = typename LinearFEM<DataTypes, ElementType>::ElementStiffness;
    using DeformationGradient = sofa::type::Mat<spatial_dimensions, spatial_dimensions, Real>;
    using RotationMatrix = sofa::type::Mat<spatial_dimensions, spatial_dimensions, Real>;

    /// the concatenation of the displacement of the 4 nodes in a single vector
    using ElementDisplacement = sofa::type::Vec<NumberOfDofsInElement, Real>;

    sofa::type::vector<RotationMatrix> m_rotations;

protected:

    void computeElementRotation(
        const std::array<Coord, NumberOfNodesInElement>& nodesPosition,
        const std::array<Coord, NumberOfNodesInElement>& nodesRestPosition,
        RotationMatrix& rotationMatrix);

public:
    using LinearFEM<DataTypes, ElementType>::LinearFEM;
    using LinearFEM<DataTypes, ElementType>::m_topology;

    static Coord computeCentroid(const std::array<Coord, NumberOfNodesInElement>& nodes);

    void addForce(VecDeriv& force, const VecCoord& position, const VecCoord& restPosition) override;
    void addDForce(VecDeriv& df, const VecDeriv& dx, Real kFactor) const override;
    void buildStiffnessMatrix(sofa::core::behavior::StiffnessMatrix::Derivative& dfdx) const override;

    void computeVonMisesStress(VonMisesStressContainer<Real>& vonMisesStressContainer,
                           const VecCoord& position,
                           const VecCoord& restPosition) const override;

    Coord translation(const std::array<Coord, NumberOfNodesInElement>& nodes) const;
};

}  // namespace elasticity

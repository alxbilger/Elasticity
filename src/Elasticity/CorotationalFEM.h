#pragma once

#include <Elasticity/LinearFEM.h>

namespace elasticity
{

template <class DataTypes, class ElementType>
class CorotationalFEM : public LinearFEM<DataTypes, ElementType>
{
    using Real = sofa::Real_t<DataTypes>;
    using VecCoord = sofa::VecCoord_t<DataTypes>;
    using Coord = sofa::Coord_t<DataTypes>;
    using TopologyElement = typename LinearFEM<DataTypes, ElementType>::TopologyElement;
    using FiniteElement = typename LinearFEM<DataTypes, ElementType>::FiniteElement;

    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;
    static constexpr sofa::Size NumberOfNodesInElement = ElementType::NumberOfNodes;
    static constexpr sofa::Size NumberOfDofsInElement = NumberOfNodesInElement * spatial_dimensions;
    static constexpr sofa::Size ElementDimension = FiniteElement::ElementDimension;

    using ElementStiffness = typename LinearFEM<DataTypes, ElementType>::ElementStiffness;
    using DeformationGradient = sofa::type::Mat<spatial_dimensions, spatial_dimensions, Real>;

    sofa::type::vector<ElementStiffness> m_rotatedStiffness;
    sofa::type::vector<sofa::type::Mat<spatial_dimensions, ElementDimension, Real>> m_restJacobians;
    sofa::type::vector<sofa::type::Quat<Real>> m_rotations;

protected:

    /**
     * Update the element stiffness matrices using the element rotation
     */
    void updateStiffnessMatrices(const VecCoord& positions, const VecCoord& restPositions) override;

    /**
     * Return the rotated element stiffness matrices
     */
    const sofa::type::vector<ElementStiffness>& stiffnessMatrices() const override;

    static void applyRotation(ElementStiffness& K, const sofa::type::Quat<Real>& rotation);

    void computeRestJacobians(const VecCoord& restPosition);

public:
    using LinearFEM<DataTypes, ElementType>::LinearFEM;
    void precomputeElementStiffness(const VecCoord& restPosition, Real youngModulus, Real poissonRatio) override;


    /**
     * Return the reference element centroid
     */
    static Coord referenceElementCentroid();

    static Coord computeCentroid(const std::array<Coord, NumberOfNodesInElement>& nodes);

    /**
     * Compute the deformation gradient at the centroid of an element
     */
    static DeformationGradient deformationGradient(const sofa::type::Mat<spatial_dimensions, NumberOfNodesInElement, Real>& nodesMatrix, const sofa::type::Mat<spatial_dimensions, ElementDimension, Real>& inverseJacobian);


};

}  // namespace elasticity

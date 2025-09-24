#pragma once

#include <Elasticity/impl/FEM/BaseFEM.h>
#include <Elasticity/finiteelement/FiniteElement[all].h>

namespace elasticity
{

template <class DataTypes>
class BaseNonLinearFEM : public BaseFEM<DataTypes>
{
};

template <class DataTypes, class ElementType>
class NonLinearFEM : public BaseNonLinearFEM<DataTypes>
{
protected:
    using VecCoord = sofa::VecCoord_t<DataTypes>;
    using VecDeriv = sofa::VecDeriv_t<DataTypes>;
    using Coord = sofa::Coord_t<DataTypes>;
    using Deriv = sofa::Deriv_t<DataTypes>;
    using Real = sofa::Real_t<DataTypes>;
    using TopologyElement = sofa::topology::Element<ElementType>;
    using FiniteElement = elasticity::FiniteElement<ElementType, DataTypes>;

    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;
    static constexpr sofa::Size NumberOfNodesInElement = ElementType::NumberOfNodes;
    static constexpr sofa::Size NumberOfDofsInElement = NumberOfNodesInElement * spatial_dimensions;
    static constexpr sofa::Size ElementDimension = FiniteElement::ElementDimension;

    using DeformationGradient = sofa::type::Mat<spatial_dimensions, spatial_dimensions, Real>;

public:
    explicit NonLinearFEM(sofa::core::topology::BaseMeshTopology* topology = nullptr);

    sofa::core::topology::BaseMeshTopology* m_topology { nullptr };

    void addForce(VecDeriv& force, const VecCoord& position, const VecCoord& restPosition) override;
    void addDForce(VecDeriv& df, const VecDeriv& dx, Real kFactor) const override;
    void buildStiffnessMatrix(sofa::core::behavior::StiffnessMatrix::Derivative& dfdx) const override;
};

}  // namespace elasticity

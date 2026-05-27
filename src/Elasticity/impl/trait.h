#pragma once

#include <sofa/fem/FiniteElement.h>
#include <sofa/component/solidmechanics/fem/elastic/impl/ElementStiffnessMatrix.h>
#include <Elasticity/impl/FullySymmetric4Tensor.h>
#include <sofa/core/trait/DataTypes.h>

namespace elasticity
{

template <class DataTypes, class ElementType>
struct trait
{
    using DataVecCoord = sofa::DataVecDeriv_t<DataTypes>;
    using DataVecDeriv = sofa::DataVecDeriv_t<DataTypes>;
    using VecCoord = sofa::VecCoord_t<DataTypes>;
    using VecDeriv = sofa::VecDeriv_t<DataTypes>;
    using Coord = sofa::Coord_t<DataTypes>;
    using Deriv = sofa::Deriv_t<DataTypes>;
    using Real = sofa::Real_t<DataTypes>;

    using FiniteElement = sofa::fem::FiniteElement<ElementType, DataTypes>;
    using TopologyElement = typename FiniteElement::TopologyElement;
    using ReferenceCoord = typename FiniteElement::ReferenceCoord;

    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;
    static constexpr sofa::Size NumberOfNodesInElement = ElementType::NumberOfNodes;
    static constexpr sofa::Size NumberOfDofsInElement = NumberOfNodesInElement * spatial_dimensions;
    static constexpr sofa::Size TopologicalDimension = FiniteElement::TopologicalDimension;
    static constexpr sofa::Size NbQuadraturePoints = FiniteElement::quadraturePoints().size();

    /// type of 2nd-order tensor for the elasticity tensor for isotropic materials
    using ElasticityTensor = FullySymmetric4Tensor<DataTypes>;

    /// the type of B in e = B d, if e is the strain, and d is the displacement
    using StrainDisplacement = sofa::component::solidmechanics::fem::elastic::StrainDisplacement<DataTypes, ElementType>;

    /// the concatenation of the displacement of the element nodes in a single vector
    using ElementDisplacement = sofa::type::Vec<NumberOfDofsInElement, Real>;

    /// tells how to compute the matrix-vector product of the stiffness matrix with a displacement
    /// vector. It does not change the result, but it can have an impact on performances.
    static constexpr sofa::component::solidmechanics::fem::elastic::MatrixVectorProductType matrixVectorProductType =
        NbQuadraturePoints > 1
            ? sofa::component::solidmechanics::fem::elastic::MatrixVectorProductType::Dense
            : sofa::component::solidmechanics::fem::elastic::MatrixVectorProductType::Factorization;

    /// the type of the element stiffness matrix
    using ElementStiffness =
        sofa::component::solidmechanics::fem::elastic::FactorizedElementStiffness<DataTypes, ElementType, matrixVectorProductType>;

    using ElementForce = sofa::type::Vec<trait::NumberOfDofsInElement, sofa::Real_t<DataTypes>>;
};

}

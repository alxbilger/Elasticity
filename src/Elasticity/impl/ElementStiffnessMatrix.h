#pragma once

#include <sofa/type/Mat.h>

#include <Elasticity/impl/ElasticityTensor.h>
#include <Elasticity/impl/MatrixTools.h>
#include <Elasticity/impl/StrainDisplacement.h>

namespace elasticity
{


template <class DataTypes, class ElementType>
struct ElementStiffnessSparse
{
    using Real = sofa::Real_t<DataTypes>;

    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;
    static constexpr sofa::Size NumberOfNodesInElement = ElementType::NumberOfNodes;
    static constexpr sofa::Size NumberOfDofsInElement = NumberOfNodesInElement * spatial_dimensions;
    static constexpr sofa::Size NumberQuadraturePoints = std::tuple_size_v<decltype(FiniteElement<ElementType, DataTypes>::quadraturePoints())>;

    static constexpr sofa::Size nbLines = NumberOfDofsInElement;
    static constexpr sofa::Size nbCols  = NumberOfDofsInElement;

    ElasticityTensorSparse<DataTypes> m_elasticityTensor;
    std::array<StrainDisplacementSparse<DataTypes, ElementType>, NumberQuadraturePoints> m_strainDisplacement;
    std::array<Real, NumberQuadraturePoints> m_factors;

    void setElasticityTensor(const ElasticityTensorSparse<DataTypes>& elasticityTensor)
    {
        m_elasticityTensor = elasticityTensor;
    }

    void setStrainDisplacement(sofa::Index quadraturePointIndex, const StrainDisplacementSparse<DataTypes, ElementType>& strainDisplacement)
    {
        m_strainDisplacement[quadraturePointIndex] = strainDisplacement;
    }

    void setFactor(sofa::Index quadraturePointIndex, Real factor)
    {
        m_factors[quadraturePointIndex] = factor;
    }

    sofa::type::Vec<NumberOfDofsInElement, Real> operator*(const sofa::type::Vec<NumberOfDofsInElement, Real>& v) const
    {
        sofa::type::Vec<NumberOfDofsInElement, Real> res;

        for (sofa::Size i = 0; i < NumberQuadraturePoints; ++i)
        {
            res += m_factors[i] * (m_strainDisplacement[i].multTranspose(m_elasticityTensor * (m_strainDisplacement[i] * v)));
        }

        return res;
    }

    void getsub(sofa::Index i, sofa::Index j, sofa::type::Mat<spatial_dimensions, spatial_dimensions, Real>& sub) const
    {

    }
};


template <class DataTypes, class ElementType>
using ElementStiffnessDense = sofa::type::Mat<
    ElementType::NumberOfNodes * DataTypes::spatial_dimensions,
    ElementType::NumberOfNodes * DataTypes::spatial_dimensions,
    sofa::Real_t<DataTypes>
>;

template <class DataTypes, class ElementType, ComputationStrategy strategy>
using ElementStiffnessMatrix = std::conditional_t<strategy == ComputationStrategy::SPARSE,
    ElementStiffnessSparse<DataTypes, ElementType>,
    ElementStiffnessDense<DataTypes, ElementType>>;


template <class DataTypes, class ElementType>
void addToStiffnessMatrix(
    sofa::Index quadraturePointIndex,
    sofa::type::Mat<
        ElementType::NumberOfNodes * DataTypes::spatial_dimensions,
        ElementType::NumberOfNodes * DataTypes::spatial_dimensions,
        sofa::Real_t<DataTypes> >& stiffnessMatrix,
    const sofa::type::Mat<
        symmetric_tensor::NumberOfIndependentElements<DataTypes::spatial_dimensions>,
        ElementType::NumberOfNodes * DataTypes::spatial_dimensions,
        sofa::Real_t<DataTypes>>& strainDisplacement,
    const sofa::type::Mat<
        symmetric_tensor::NumberOfIndependentElements<DataTypes::spatial_dimensions>,
        symmetric_tensor::NumberOfIndependentElements<DataTypes::spatial_dimensions>,
        sofa::Real_t<DataTypes>>& elasticityTensor,
    sofa::Real_t<DataTypes> factor)
{
    stiffnessMatrix += factor * strainDisplacement.transposed() * elasticityTensor * strainDisplacement;
}

template <class DataTypes, class ElementType>
void addToStiffnessMatrix(
    sofa::Index quadraturePointIndex,
    ElementStiffnessSparse<DataTypes, ElementType>& stiffness,
    const StrainDisplacementSparse<DataTypes, ElementType>& strainDisplacement,
    const ElasticityTensorSparse<DataTypes>& elasticityTensor,
    sofa::Real_t<DataTypes> factor)
{
    stiffness.setStrainDisplacement(quadraturePointIndex, strainDisplacement);
    stiffness.setElasticityTensor(elasticityTensor);
    stiffness.setFactor(quadraturePointIndex, factor);
}


template <class DataTypes, class ElementType, ComputationStrategy strategy>
ElementStiffnessMatrix<DataTypes, ElementType, strategy> integrate(
    const std::array<sofa::Coord_t<DataTypes>, ElementType::NumberOfNodes>& nodesCoordinates,
    const ElasticityTensor<DataTypes, strategy>& elasticityTensor)
{
    using Real = sofa::Real_t<DataTypes>;
    using FiniteElement = FiniteElement<ElementType, DataTypes>;

    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;
    static constexpr sofa::Size NumberOfNodesInElement = ElementType::NumberOfNodes;
    static constexpr sofa::Size ElementDimension = FiniteElement::ElementDimension;

    ElementStiffnessMatrix<DataTypes, ElementType, strategy> K;

    sofa::Size quadraturePointIndex = 0;
    for (const auto& [quadraturePoint, weight] : FiniteElement::quadraturePoints())
    {
        // gradient of shape functions in the reference element evaluated at the quadrature point
        const sofa::type::Mat<NumberOfNodesInElement, ElementDimension, Real> dN_dq_ref =
            FiniteElement::gradientShapeFunctions(quadraturePoint);

        // jacobian of the mapping from the reference space to the physical space, evaluated at the
        // quadrature point
        sofa::type::Mat<spatial_dimensions, ElementDimension, Real> jacobian;
        for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
            jacobian += sofa::type::dyad(nodesCoordinates[i], dN_dq_ref[i]);

        const auto detJ = elasticity::determinant(jacobian);
        const sofa::type::Mat<ElementDimension, spatial_dimensions, Real> J_inv =
            elasticity::inverse(jacobian);

        // gradient of the shape functions in the physical element evaluated at the quadrature point
        sofa::type::Mat<NumberOfNodesInElement, spatial_dimensions, Real> dN_dq(sofa::type::NOINIT);
        for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
            dN_dq[i] = J_inv.transposed() * dN_dq_ref[i];

        const auto B = makeStrainDisplacement<DataTypes, ElementType, strategy>(dN_dq);
        addToStiffnessMatrix<DataTypes, ElementType>(quadraturePointIndex++, K, B, elasticityTensor, weight * detJ);
    }
    return K;
}



}

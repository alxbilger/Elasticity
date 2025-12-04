#pragma once

#include <Elasticity/finiteelement/FiniteElement.h>
#include <Elasticity/impl/FullySymmetric4Tensor.h>
#include <Elasticity/impl/MatrixTools.h>
#include <Elasticity/impl/StrainDisplacement.h>
#include <sofa/helper/MatEigen.h>
#include <sofa/type/Mat.h>

namespace elasticity
{

template <class DataTypes, class ElementType>
using ElementStiffness = sofa::type::MatSym<
    ElementType::NumberOfNodes * DataTypes::spatial_dimensions,
    sofa::Real_t<DataTypes>
>;

enum class MatrixVectorProductType
{
    Factorization,
    Dense
};

template <class DataTypes, class ElementType, MatrixVectorProductType matrixVectorProductType>
struct FactorizedElementStiffness
{
private:
    using Real = sofa::Real_t<DataTypes>;
    using FiniteElement = elasticity::FiniteElement<ElementType, DataTypes>;
    static constexpr auto NbQuadraturePoints = FiniteElement::quadraturePoints().size();


    /// The factors of the stiffness matrix
    /// @{
    std::array<StrainDisplacement<DataTypes, ElementType>, NbQuadraturePoints> B;
    FullySymmetric4Tensor<DataTypes> elasticityTensor;
    static constexpr sofa::Size NumberOfIndependentElements = symmetric_tensor::NumberOfIndependentElements<DataTypes::spatial_dimensions>;
    sofa::type::Mat<NumberOfIndependentElements, NumberOfIndependentElements, Real> elasticityTensorMat;
    std::array<Real, NbQuadraturePoints> factors;
    /// @}

    /**
     * The assembled stiffness matrix
     */
    sofa::type::Mat<
        ElementType::NumberOfNodes * DataTypes::spatial_dimensions,
        ElementType::NumberOfNodes * DataTypes::spatial_dimensions,
        sofa::Real_t<DataTypes>> stiffnessMatrix;

    using Vec = sofa::type::Vec<
        ElementType::NumberOfNodes * DataTypes::spatial_dimensions,
        sofa::Real_t<DataTypes>>;

public:
    void setElasticityTensor(const FullySymmetric4Tensor<DataTypes>& elasticityTensor)
    {
        this->elasticityTensor = elasticityTensor;
        for (std::size_t i = 0; i < NumberOfIndependentElements; ++i)
        {
            for (std::size_t j = 0; j < NumberOfIndependentElements; ++j)
            {
                elasticityTensorMat(i, j) = elasticityTensor.toVoigtMatSym()(i, j);
            }
        }
    }

    void add(std::size_t quadraturePointIndex,
        const Real factor,
        const StrainDisplacement<DataTypes, ElementType>& B_)
    {
        this->factors[quadraturePointIndex] = factor;
        this->B[quadraturePointIndex] = B_;

        const auto K_i = factor * B_.multTranspose(elasticityTensorMat * B_);
        stiffnessMatrix += K_i;
    }


    Vec operator*(const Vec& v) const
    {
        if constexpr (matrixVectorProductType == MatrixVectorProductType::Factorization)
        {
            Vec result;
            for (std::size_t i = 0; i < NbQuadraturePoints; ++i)
            {
                const auto Bv = B[i] * v;
                const auto CBv = elasticityTensorMat * Bv;
                const auto BTCBv = B[i].multTranspose(CBv);
                result += factors[i] * BTCBv;
            }
            return result;
        }
        else
        {
            return stiffnessMatrix * v;
        }
    }

    const auto& getAssembledMatrix() const { return stiffnessMatrix; }
};

template <class DataTypes, class ElementType, MatrixVectorProductType matrixVectorProductType = MatrixVectorProductType::Dense>
FactorizedElementStiffness<DataTypes, ElementType, matrixVectorProductType> integrate(
    const std::array<sofa::Coord_t<DataTypes>, ElementType::NumberOfNodes>& nodesCoordinates,
    const FullySymmetric4Tensor<DataTypes>& elasticityTensor)
{
    using Real = sofa::Real_t<DataTypes>;
    using FiniteElement = FiniteElement<ElementType, DataTypes>;

    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;
    static constexpr sofa::Size NumberOfNodesInElement = ElementType::NumberOfNodes;
    static constexpr sofa::Size TopologicalDimension = FiniteElement::TopologicalDimension;

    FactorizedElementStiffness<DataTypes, ElementType, matrixVectorProductType> K;
    K.setElasticityTensor(elasticityTensor);

    std::size_t quadraturePointIndex = 0;
    for (const auto& [quadraturePoint, weight] : FiniteElement::quadraturePoints())
    {
        // gradient of shape functions in the reference element evaluated at the quadrature point
        const sofa::type::Mat<NumberOfNodesInElement, TopologicalDimension, Real> dN_dq_ref =
            FiniteElement::gradientShapeFunctions(quadraturePoint);

        // jacobian of the mapping from the reference space to the physical space, evaluated at the
        // quadrature point
        sofa::type::Mat<spatial_dimensions, TopologicalDimension, Real> jacobian;
        for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
            jacobian += sofa::type::dyad(nodesCoordinates[i], dN_dq_ref[i]);

        const auto detJ = elasticity::determinant(jacobian);
        const sofa::type::Mat<TopologicalDimension, spatial_dimensions, Real> J_inv =
            elasticity::inverse(jacobian);

        // gradient of the shape functions in the physical element evaluated at the quadrature point
        sofa::type::Mat<NumberOfNodesInElement, spatial_dimensions, Real> dN_dq(sofa::type::NOINIT);
        for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
            dN_dq[i] = J_inv.transposed() * dN_dq_ref[i];

        const auto B = makeStrainDisplacement<DataTypes, ElementType>(dN_dq);
        K.add(quadraturePointIndex, weight * detJ, B);

        ++quadraturePointIndex;
    }
    return K;
}



}

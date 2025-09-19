#pragma once

#include <Elasticity/impl/SymmetricTensor.h>
#include <sofa/core/trait/DataTypes.h>
#include <sofa/type/Mat.h>

namespace elasticity
{

/**
 * @brief Converts Young's modulus and Poisson's ratio to Lamé parameters.
 *
 * This function calculates and returns the two Lamé parameters, μ (shear modulus)
 * and λ, derived from the given Young’s modulus and Poisson’s ratio of a material.
 * These parameters are fundamental in describing isotropic elastic behavior.
 *
 * @param youngModulus The Young's modulus of the material, representing its stiffness.
 * @param poissonRatio The Poisson's ratio of the material, describing its deformation behavior.
 * @return A pair containing the calculated Lamé parameters:
 *         - First: μ (shear modulus),
 *         - Second: λ.
 */
template<class DataTypes>
std::pair<sofa::Real_t<DataTypes>, sofa::Real_t<DataTypes>>
toLameParameters(sofa::Real_t<DataTypes> youngModulus, sofa::Real_t<DataTypes> poissonRatio)
{
    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;
    const auto mu = youngModulus / (2 * (1 + poissonRatio));
    const auto lambda = youngModulus * poissonRatio / ((1 + poissonRatio) * (1 - (spatial_dimensions - 1) * poissonRatio));
    return std::make_pair(mu, lambda);
}

template <class DataTypes>
using ElasticityTensor = sofa::type::Mat<
    symmetric_tensor::NumberOfIndependentElements<DataTypes::spatial_dimensions>,
    symmetric_tensor::NumberOfIndependentElements<DataTypes::spatial_dimensions>,
    sofa::Real_t<DataTypes>>;

/**
 * @brief Creates an isotropic elasticity tensor for given material properties.
 *
 * This function constructs and returns an elasticity tensor for an isotropic material
 * characterized by its Young's modulus and Poisson's ratio. It computes the tensor
 * using the Lamé parameters, which are derived from the given material properties.
 *
 * @param youngModulus The Young's modulus of the material, representing its stiffness.
 * @param poissonRatio The Poisson's ratio of the material, describing its deformation behavior.
 * @return The isotropic elasticity tensor represented as an object of
 * `ElasticityTensor<DataTypes>`.
 */
template <class DataTypes>
ElasticityTensor<DataTypes> makeIsotropicElasticityTensor(sofa::Real_t<DataTypes> youngModulus, sofa::Real_t<DataTypes> poissonRatio)
{
    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;
    static constexpr sofa::Size NumberOfIndependentElements = symmetric_tensor::NumberOfIndependentElements<spatial_dimensions>;
    using Real = sofa::Real_t<DataTypes>;

    const auto [mu, lambda] = toLameParameters<DataTypes>(youngModulus, poissonRatio);
    ElasticityTensor<DataTypes> C;

    for (sofa::Size i = 0; i < spatial_dimensions; ++i)
    {
        for (sofa::Size j = 0; j < spatial_dimensions; ++j)
        {
            C(i, j) = lambda;
        }
        C(i, i) += 2 * mu;
    }

    for (sofa::Size i = spatial_dimensions; i < NumberOfIndependentElements; ++i)
    {
        C(i, i) = mu;
    }

    return C;
}

}

#pragma once

#include <Elasticity/impl/SymmetricTensor.h>
#include <sofa/core/trait/DataTypes.h>
#include <sofa/type/Mat.h>
#include <sofa/type/MatSym.h>

#include <Elasticity/impl/KroneckerDelta.h>
#include <Elasticity/impl/VoigtNotation.h>

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

/**
 * A class to represent the Lagrangian elasticity tensor.
 *
 * The elasticity tensor is the derivative of the second Piola-Kirchhoff with respect to the
 * Green-Lagrange tensor. It is a 4th-order tensor
 */
template <class DataTypes>
class ElasticityTensor
{
private:
    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;
    static constexpr sofa::Size NumberOfIndependentElements = symmetric_tensor::NumberOfIndependentElements<spatial_dimensions>;
    using Real = sofa::Real_t<DataTypes>;

public:
    ElasticityTensor() = default;

    template<class Callable>
    explicit ElasticityTensor(Callable callable) : m_matrix(sofa::type::NOINIT)
    {
        fill(callable);
    }

    template<class Callable>
    void fill(Callable callable)
    {
#ifndef NDEBUG
        checkSymmetry(callable);
#endif

        for (sofa::Size a = 0; a < NumberOfIndependentElements; ++a)
        {
            const auto [i, j] = voigtIndices<DataTypes>(a);
            for (sofa::Size b = a; b < NumberOfIndependentElements; ++b) // the Voigt representation is symmetric, that is why j starts at i
            {
                const auto [k, l] = voigtIndices<DataTypes>(b);
                m_matrix(a, b) = callable(i, j, k, l);
            }
        }
    }

    Real& operator()(sofa::Size i, sofa::Size j, sofa::Size k, sofa::Size l)
    {
        const auto a = voigtIndex<DataTypes>(i, j);
        const auto b = voigtIndex<DataTypes>(k, l);
        return m_matrix(a, b);
    }

    Real operator()(sofa::Size i, sofa::Size j, sofa::Size k, sofa::Size l) const
    {
        const auto a = voigtIndex<DataTypes>(i, j);
        const auto b = voigtIndex<DataTypes>(k, l);
        return m_matrix(a, b);
    }

    const sofa::type::MatSym<NumberOfIndependentElements, Real>& toVoigtMatSym() const
    {
        return m_matrix;
    }

private:
    sofa::type::MatSym<NumberOfIndependentElements, Real> m_matrix;

    template<class Callable>
    static void checkSymmetry(Callable callable)
    {
        for (sofa::Size i = 0; i < spatial_dimensions; ++i)
        {
            for (sofa::Size j = 0; j < spatial_dimensions; ++j)
            {
                for (sofa::Size k = 0; k < spatial_dimensions; ++k)
                {
                    for (sofa::Size l = 0; l < spatial_dimensions; ++l)
                    {
                        const auto ijkl = callable(i, j, k, l);
                        const auto klij = callable(k, l, i, j);
                        const auto ijlk = callable(i, j, l, k);
                        assert(std::abs(ijkl - klij) < 1e-12 );
                        assert(std::abs(ijkl - ijlk) < 1e-12 );
                    }
                }
            }
        }
    }
};


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
    ElasticityTensor<DataTypes> C(
        [mu, lambda](sofa::Index i, sofa::Index j, sofa::Index k, sofa::Index l)
        {
            return mu * (kroneckerDelta<Real>(i, k) * kroneckerDelta<Real>(j, l) + kroneckerDelta<Real>(i, l) * kroneckerDelta<Real>(j, k)) +
                        lambda * kroneckerDelta<Real>(i, j) * kroneckerDelta<Real>(k, l);
        });

    return C;
}

}

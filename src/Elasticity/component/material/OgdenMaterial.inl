#pragma once

#include <Elasticity/component/material/OgdenMaterial.h>
#include <Elasticity/component/PK2HyperelasticMaterial.inl>
#include <Eigen/Eigenvalues>

namespace elasticity
{

template <class DataTypes>
OgdenMaterial<DataTypes>::OgdenMaterial()
    : m_mu(initData(&m_mu, static_cast<Real>(1e2), "mu",
                      "Material constant relevant to the shear modulus"))
    , m_alpha(initData(&m_alpha, static_cast<Real>(1.5), "alpha",
                      "Material constant exponent related to strain-stiffening"))
    , m_kappa(initData(&m_kappa, static_cast<Real>(1e3), "kappa",
                      "Material constant related to the bulk modulus"))
{}

template <class DataTypes>
auto OgdenMaterial<DataTypes>::secondPiolaKirchhoffStress(Strain<DataTypes>& strain) -> StressTensor
{
    using EigenMatrix = Eigen::Matrix<Real, spatial_dimensions, spatial_dimensions>;

    const auto& C = strain.getRightCauchyGreenTensor();
    const auto C_1 = inverse(C);

    const auto J = strain.getDeterminantDeformationGradient();
    assert(J > 0);

    const auto S_isochoric = [this, J, &strain, &C_1, &C] ()
    {
        const Real mu = m_mu.getValue();
        const Real alpha = m_alpha.getValue();

        const Real FJ = pow(J, -alpha / static_cast<Real>(3));

        EigenMatrix CEigen;
        for (sofa::Index m = 0; m < spatial_dimensions; ++m)
            for (sofa::Index n = 0; n < spatial_dimensions; ++n) 
                CEigen(m, n) = C(m, n);

        // Disable temporarilly until fixed /*Eigen::SelfAdjointEigenSolver<EigenMatrix>*/
        Eigen::EigenSolver<EigenMatrix> EigenProblemSolver(CEigen, true);
        if (EigenProblemSolver.info() != Eigen::Success)
            dmsg_warning("OgdenMaterial") << "EigenSolver iterations failed to converge";
        const auto eigenvectors = EigenProblemSolver.eigenvectors().real();
        const auto eigenvalues = EigenProblemSolver.eigenvalues().real();

        // Precompute trace of C^(alpha/2) and eigenvalue powers: lambda^(alpha/2 - 1)
        const Real aBy2 = static_cast<Real>(0.5) * alpha;
        const Real aBy2Minus1 = aBy2 - static_cast<Real>(1);
        Real trCaBy2{static_cast<Real>(0)};
        typename DataTypes::Coord eigenvaluePowers{};
        for (sofa::Index n = 0; n < spatial_dimensions; ++n)
        {
            trCaBy2 += pow(eigenvalues(n), aBy2);
            eigenvaluePowers(n) = pow(eigenvalues(n), aBy2Minus1);
        }

        const Real coefS1 = FJ * mu / (alpha * static_cast<Real>(2));
        const Real coefS2 = -FJ * mu / (alpha * static_cast<Real>(6)) * trCaBy2;

        StressTensor S{};

        for (sofa::Index i = 0; i < spatial_dimensions; ++i)
        {
            for (sofa::Index j = 0; j < spatial_dimensions; ++j)
            {
                for (sofa::Index n = 0; n < spatial_dimensions; ++n)
                    S(i,j) += coefS1 * eigenvaluePowers(n) * eigenvectors(i,n) * eigenvectors(j,n);
                S(i, j) += coefS2 * C_1(j, i);
            }
        }

        return static_cast<Real>(2) * S;
    }();

    const auto S_volumetric = [this, J, &C_1]()
    {
        const auto kappa = m_kappa.getValue();
        return kappa * log(J) * C_1;
    }();
    
    return S_isochoric + S_volumetric;
}

template <class DataTypes>
auto OgdenMaterial<DataTypes>::elasticityTensor(Strain<DataTypes>& strain) -> ElasticityTensor
{
    return elasticity::FullySymmetric4Tensor<sofa::defaulttype::StdVectorTypes<DataTypes::Coord,DataTypes::Coord,Real>>();
}

}  // namespace elasticity

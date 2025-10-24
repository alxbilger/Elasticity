#include <Elasticity/component/material/StVenantKirchhoffMaterial.h>
#include <gtest/gtest.h>
#include <sofa/testing/LinearCongruentialRandomGenerator.h>
#include <Eigen/Eigenvalues>

namespace elasticity
{

template <typename T>
class BaseMaterialTest : public testing::Test
{
protected:
    using DataTypes = typename T::DataTypes;
    using Real = sofa::Real_t<DataTypes>;
    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;
    using DeformationGradient = sofa::type::Mat<spatial_dimensions, spatial_dimensions, Real>;

    void SetUp() override
    {
        material = sofa::core::objectmodel::New<T>();
        material->init();
    }

    DeformationGradient generatePositiveDefiniteMatrix()
    {
        DeformationGradient F;

        for (int j = 0; j < spatial_dimensions; j++)
        {
            for (int i = 0; i < spatial_dimensions; i++)
            {
                F(i, j) = lcg.generateInRange(-10., 10.);
            }
        }

        F = 0.5 * (F + F.transposed()); // Ensuring symmetry

        Eigen::Map<Eigen::Matrix<Real, spatial_dimensions, spatial_dimensions, Eigen::RowMajor>> Fmap(const_cast<Real*>(F.data()));

        // Making the matrix positive definite
        const Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Real, spatial_dimensions, spatial_dimensions> > eigensolver(Fmap);
        const Eigen::Matrix<Real, spatial_dimensions, 1>& eigenvalues = eigensolver.eigenvalues();
        const Real min_eigval = eigenvalues.minCoeff();
        if (min_eigval <= 0)
        {
            Fmap += Eigen::Matrix<Real,spatial_dimensions, spatial_dimensions>::Identity() * (-min_eigval + 1e-6); // Adding a small value to ensure positive definiteness
        }

        return F;
    }


    T::SPtr material;
    sofa::testing::LinearCongruentialRandomGenerator lcg { 31321 };
};

template <typename T>
class MaterialTest : public BaseMaterialTest<T>
{
    using BaseMaterialTest<T>::material;

public:

    using DataTypes = typename T::DataTypes;
    using Real = sofa::Real_t<DataTypes>;
    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;
    using DeformationGradient = sofa::type::Mat<spatial_dimensions, spatial_dimensions, Real>;

    void testDerivativePK1()
    {
        const auto F = this->generatePositiveDefiniteMatrix();

        const auto P = material->firstPiolaKirchhoffStress(F);
        const auto A = material->materialTangentModulus(F);

        EXPECT_TRUE(std::none_of(P.data(), P.data() + spatial_dimensions * spatial_dimensions, [](Real x) { return x == 0; })) << P;

        // Small perturbation for finite difference
        constexpr Real epsilon = 1e-12;

        std::stringstream ss;
        ss << "numerical,analytical\n";

        // For each component of F
        for(sofa::Size i = 0; i < spatial_dimensions; ++i)
        {
            for(sofa::Size j = 0; j < spatial_dimensions; ++j)
            {
                // Create perturbed deformation gradient
                auto F_perturbed = F;
                F_perturbed(i,j) += epsilon;

                // Compute perturbed stress
                const auto P_perturbed = material->firstPiolaKirchhoffStress(F_perturbed);

                // Compute numerical derivative
                const auto dP = (P_perturbed - P) / (epsilon);

                // For each component of P
                for(sofa::Size k = 0; k < spatial_dimensions; ++k)
                {
                    for(sofa::Size l = 0; l < spatial_dimensions; ++l)
                    {
                        // Compare numerical and analytical derivatives
                        const Real numerical = dP(k,l);
                        const Real analytical = A(k, l, i, j);
                        ss << numerical << "," << analytical << std::endl;
                        EXPECT_NEAR(numerical, analytical, 1e-3);
                    }
                }
            }
        }

        std::cout << ss.str();
    }

    void testMajorSymmetry()
    {
        const auto F = this->generatePositiveDefiniteMatrix();
        const auto A = material->materialTangentModulus(F);

        for(sofa::Size i = 0; i < spatial_dimensions; ++i)
        {
            for(sofa::Size j = 0; j < spatial_dimensions; ++j)
            {
                for(sofa::Size k = 0; k < spatial_dimensions; ++k)
                {
                    for(sofa::Size l = 0; l < spatial_dimensions; ++l)
                    {
                        EXPECT_NEAR(A(i, j, k, l), A(k, l, i, j), 1e-6);
                    }
                }
            }
        }
    }
};


using AllMaterials = ::testing::Types<
    StVenantKirchhoffMaterial<sofa::defaulttype::Vec3Types>,
    StVenantKirchhoffMaterial<sofa::defaulttype::Vec2Types>,
    StVenantKirchhoffMaterial<sofa::defaulttype::Vec1Types>
>;
TYPED_TEST_SUITE(MaterialTest, AllMaterials);




TYPED_TEST(MaterialTest, derivativePK1)
{
    this->testDerivativePK1();
}

TYPED_TEST(MaterialTest, tensorMajorSymmetry)
{
    this->testMajorSymmetry();
}

template <typename T>
class PK2Material : public T
{
public:
    SOFA_CLASS(PK2Material, T);
    using T::secondPiolaKirchhoffStress;
    using T::elasticityTensor;
};

template <typename T>
class PK2MaterialTest : public BaseMaterialTest<PK2Material<T>>
{
public:
    using DataTypes = typename T::DataTypes;
    using Real = sofa::Real_t<DataTypes>;
    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;

    void testMinorSymmetryPK2()
    {
        const auto F = this->generatePositiveDefiniteMatrix();
        const auto S = this->material->secondPiolaKirchhoffStress(F.transposed() * F);

        for(sofa::Size i = 0; i < spatial_dimensions; ++i)
        {
            for(sofa::Size j = 0; j < spatial_dimensions; ++j)
            {
                EXPECT_NEAR(S(i,j), S(j,i), 1e-6);
            }
        }

    }

    void testSymmetryElasticityTensor()
    {
        const auto F = this->generatePositiveDefiniteMatrix();
        const auto C = this->material->elasticityTensor(F);

        for(sofa::Size i = 0; i < spatial_dimensions; ++i)
        {
            for(sofa::Size j = 0; j < spatial_dimensions; ++j)
            {
                for(sofa::Size k = 0; k < spatial_dimensions; ++k)
                {
                    for(sofa::Size l = 0; l < spatial_dimensions; ++l)
                    {
                        EXPECT_NEAR(C(i, j, k, l), C(k, l, i, j), 1e-6);
                        EXPECT_NEAR(C(i, j, k, l), C(j, i, k, l), 1e-6);
                        EXPECT_NEAR(C(i, j, k, l), C(i, j, l, k), 1e-6);
                    }
                }
            }
        }
    }

    void testDerivativePK2()
    {
        const auto FTF = [this]()
        {
            const auto F = this->generatePositiveDefiniteMatrix();
            return F.transposed() * F;
        }();

        const auto S = this->material->secondPiolaKirchhoffStress(FTF);
        const auto C = this->material->elasticityTensor(FTF);

        // Small perturbation for finite difference
        constexpr Real epsilon = 1e-12;

        std::stringstream ss;
        ss << "numerical,analytical\n";

        // For each component of F
        for(sofa::Size i = 0; i < spatial_dimensions; ++i)
        {
            for(sofa::Size j = 0; j < spatial_dimensions; ++j)
            {
                // Create perturbed deformation gradient
                auto FTF_perturbed = FTF;

                // For symmetric tensor, perturb both (i,j) and (j,i)
                if (i == j)
                {
                    FTF_perturbed(i,j) += epsilon;
                }
                else
                {
                    FTF_perturbed(i,j) += epsilon;
                    FTF_perturbed(j,i) += epsilon;
                }

                // Compute perturbed stress
                const auto S_perturbed = this->material->secondPiolaKirchhoffStress(FTF_perturbed);

                // Compute numerical derivative
                // For off-diagonal, we perturbed by 2*epsilon total (both (i,j) and (j,i))
                const auto dS = (i == j) ? (S_perturbed - S) / epsilon : (S_perturbed - S) / (2 * epsilon);

                // For each component of S
                for(sofa::Size k = 0; k < spatial_dimensions; ++k)
                {
                    for(sofa::Size l = 0; l < spatial_dimensions; ++l)
                    {
                        // Compare numerical and analytical derivatives
                        const Real numerical = dS(k,l);
                        // For off-diagonal terms, the analytical derivative needs to account for symmetry
                        const Real analytical = (i == j) ? C(i, j, k, l) : (C(i, j, k, l) + C(j, i, k, l));
                        ss << numerical << "," << analytical << std::endl;
                        EXPECT_NEAR(numerical, analytical, 1e-3);
                    }
                }
            }
        }

        std::cout << ss.str();
    }
    
};

using PK2Materials = ::testing::Types<
    StVenantKirchhoffMaterial<sofa::defaulttype::Vec3Types>,
    StVenantKirchhoffMaterial<sofa::defaulttype::Vec2Types>,
    StVenantKirchhoffMaterial<sofa::defaulttype::Vec1Types>
>;
TYPED_TEST_SUITE(PK2MaterialTest, PK2Materials);

TYPED_TEST(PK2MaterialTest, symmetryElasticityTensor)
{
    this->testSymmetryElasticityTensor();
}

TYPED_TEST(PK2MaterialTest, minorSymmetryPK2)
{
    this->testMinorSymmetryPK2();
}

TYPED_TEST(PK2MaterialTest, derivativePK2)
{
    this->testDerivativePK2();
}

}

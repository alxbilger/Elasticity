#include <Elasticity/component/material/StVenantKirchhoffMaterial.h>
#include <gtest/gtest.h>
#include <sofa/testing/LinearCongruentialRandomGenerator.h>
#include <Eigen/Eigenvalues>

namespace elasticity
{

template <typename T>
class MaterialTest : public testing::Test
{
public:

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

    void testDerivative()
    {
        const auto F = generatePositiveDefiniteMatrix();

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


    T::SPtr material;
    sofa::testing::LinearCongruentialRandomGenerator lcg { 31321 };
};


using AllMaterials = ::testing::Types<
    StVenantKirchhoffMaterial<sofa::defaulttype::Vec3Types>,
    StVenantKirchhoffMaterial<sofa::defaulttype::Vec2Types>,
    StVenantKirchhoffMaterial<sofa::defaulttype::Vec1Types>
>;
TYPED_TEST_SUITE(MaterialTest, AllMaterials);




TYPED_TEST(MaterialTest, derivativeConsistency)
{
    this->testDerivative();
}


}

#pragma once

#include <Elasticity/component/material/YeohThirdOrderMaterial.h>
#include <Elasticity/component/PK2HyperelasticMaterial.inl>

namespace elasticity
{

template <class DataTypes>
YeohThirdOrderMaterial<DataTypes>::YeohThirdOrderMaterial()
    : m_c1(initData(&m_c1, static_cast<Real>(1e3), "C1",
                      "First material constant"))
    , m_c2(initData(&m_c2, static_cast<Real>(1e3), "C2",
                      "Second material constant"))
    , m_c3(initData(&m_c3, static_cast<Real>(1e3), "C3",
                      "Third material constant"))
    , m_bulkModulus1(initData(&m_bulkModulus1, static_cast<Real>(1e2), "bulkModulus1", "First bulk modulus"))
    , m_bulkModulus2(initData(&m_bulkModulus2, static_cast<Real>(1e2), "bulkModulus2", "Second bulk modulus"))
    , m_bulkModulus3(initData(&m_bulkModulus3, static_cast<Real>(1e2), "bulkModulus3", "Third bulk modulus"))
{}

template <class DataTypes>
auto YeohThirdOrderMaterial<DataTypes>::secondPiolaKirchhoffStress(Strain<DataTypes>& strain) -> StressTensor
{
    const auto& C = strain.getRightCauchyGreenTensor();
    const auto C_1 = inverse(C);

    const auto J = strain.getDeterminantDeformationGradient();
    assert(J > 0);

    const auto S_isochoric = [this, J, &strain, &C_1, &C]()
    {
        static constexpr auto& I = Strain<DataTypes>::identity;
        static constexpr Real dim_1 = static_cast<Real>(1) / static_cast<Real>(spatial_dimensions);

        const auto c1 = m_c1.getValue();
        const auto c2 = m_c2.getValue();
        const auto c3 = m_c3.getValue();

        const auto invariant1 = strain.getInvariant1();

        const auto common_part = pow(J,-static_cast<Real>(2) * dim_1)*(I-dim_1 *invariant1*C_1);

        const auto S_c1 = c1;
        const auto S_c2 = 2*c2*(pow(J,-static_cast<Real>(2) * dim_1)*invariant1-3);
        const auto S_c3 = 3*c3*pow((pow(J,-static_cast<Real>(2) * dim_1)*invariant1-3),2);

        return static_cast<Real>(2) * common_part* (S_c1+S_c2+S_c3);
    }();

    const auto S_volumetric = [this, J, &C_1]()
    {
        const auto bulk1 = m_bulkModulus1.getValue();
        const auto bulk2 = m_bulkModulus2.getValue();
        const auto bulk3 = m_bulkModulus3.getValue();

        return  C_1*(bulk1*log(J)+2*bulk2*pow(log(J),3)+3*bulk3*pow(log(J),5));
    }();

    return S_isochoric + S_volumetric;
}

template <class DataTypes>
auto YeohThirdOrderMaterial<DataTypes>::elasticityTensor(Strain<DataTypes>& strain) -> ElasticityTensor
{
    static constexpr Real dim_1 = static_cast<Real>(1) / static_cast<Real>(spatial_dimensions);
    const auto& C = strain.getRightCauchyGreenTensor();
    const auto J = strain.getDeterminantDeformationGradient();
    const auto logJ = log(J);
    const auto J_2dim = pow(J, -2 * dim_1);
    const auto J_4dim = pow(J, -4 * dim_1);
    const auto C_1 = elasticity::inverse(C);
    const auto I1 = strain.getInvariant1();

    const auto c1 = m_c1.getValue();
    const auto c2 = m_c2.getValue();
    const auto c3 = m_c3.getValue();
    const auto bulk1 = m_bulkModulus1.getValue();
    const auto bulk2 = m_bulkModulus2.getValue();
    const auto bulk3 = m_bulkModulus3.getValue();

 return ElasticityTensor(
        [&](sofa::Index i, sofa::Index j, sofa::Index k, sofa::Index l)
        {
            //derivative of C^{-1} with respect to C
            const auto dC_1dC = -static_cast<Real>(0.5) * (C_1(i, k) * C_1(l, j) + C_1(i, l) * C_1(k, j));
            
            const Real dS_part1 = 2*(c1+2*c2*(J_2dim*I1-3)+3*c3*pow(J_2dim*I1-3,2))*(-J_2dim*dim_1*C_1(k,l)*(kroneckerDelta<Real>(i, j)-I1*dim_1*C_1(i,j))+J_2dim*dim_1*(-C_1(i,j)*kroneckerDelta<Real>(k, l)+I1*dC_1dC));
            

            const Real dS_part2 = 2*J_4dim*(kroneckerDelta<Real>(i, j)-I1*dim_1*C_1(i,j))*(kroneckerDelta<Real>(k, l)-I1*dim_1*C_1(k,l))*(2*c2+6*c3*(J_2dim*I1-3));

            const Real dS_isochoric_dC = dS_part1+dS_part2;

            // the derivative of S_volumetric with respect to C
            // this term has both minor and major symmetries
            const Real dS_volumetric_dC = dC_1dC * (bulk1 *logJ + 2 * bulk2 * pow(logJ,3)+3 * bulk3 *pow(logJ,5))  + 0.5 * C_1(l, k) * C_1(i, j)* (bulk1 + 6 * bulk2 * pow(logJ,2)+15 * bulk3 * pow(logJ,4));

            return 2 * (dS_isochoric_dC + dS_volumetric_dC);
        });
}

}  // namespace elasticity

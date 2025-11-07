#pragma once

#include <sofa/type/Mat.h>
#include <Elasticity/impl/MatrixTools.h>

namespace elasticity
{

inline struct DeformationGradientTag {} deformationGradient;
inline struct RightCauchyGreenTensorTag {} rightCauchyGreenTensor;

template <class DataTypes>
struct Strain
{
    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;
    using Real = sofa::Real_t<DataTypes>;
    using DeformationGradient = sofa::type::Mat<spatial_dimensions, spatial_dimensions, Real>;
    using RightCauchyGreenTensor = sofa::type::Mat<spatial_dimensions, spatial_dimensions, Real>;
    using GreenLagrangeTensor = sofa::type::Mat<spatial_dimensions, spatial_dimensions, Real>;
    static constexpr IdentityMatrix identity;

    constexpr Strain(DeformationGradientTag, const DeformationGradient& F) : m_deformationGradient(F) {}
    constexpr Strain(RightCauchyGreenTensorTag, const DeformationGradient& C) : m_rightCauchyGreenTensor(C) {}

    const DeformationGradient& deformationGradient() const
    {
        if (!m_deformationGradient.has_value())
        {
            throw std::logic_error("Cannot access deformation gradient. It has not been provided");
        }
        return *m_deformationGradient;
    }

    void setDeformationGradient(const DeformationGradient& F)
    {
        reset();
        m_deformationGradient = F;
    }

    Real getDeterminantDeformationGradient()
    {
        if (!m_determinantF)
        {
            computeDeterminant();
        }
        return *m_determinantF;
    }

    const RightCauchyGreenTensor& getRightCauchyGreenTensor()
    {
        if (!m_rightCauchyGreenTensor)
        {
            computeRightCauchyGreenTensor();
        }
        return *m_rightCauchyGreenTensor;
    }

    void setRightCauchyGreenTensor(const RightCauchyGreenTensor& C)
    {
        reset();
        m_rightCauchyGreenTensor = C;
    }

    const GreenLagrangeTensor& getGreenLagrangeTensor()
    {
        if (!m_greenLagrangeTensor)
        {
            computeGreenLagrangeTensor();
        }
        return *m_greenLagrangeTensor;
    }

    Real getInvariant1()
    {
        if (!m_invariant1)
        {
            computeInvariant1();
        }
        return *m_invariant1;
    }

    Real getInvariant2()
    {
        if (!m_invariant2)
        {
            computeInvariant2();
        }
        return *m_invariant2;
    }

    Real getInvariant3()
    {
        if (!m_invariant3)
        {
            computeInvariant3();
        }
        return *m_invariant3;
    }

    void reset()
    {
        m_deformationGradient.reset();
        m_rightCauchyGreenTensor.reset();
        m_determinantF.reset();
        m_greenLagrangeTensor.reset();
        m_invariant1.reset();
        m_invariant2.reset();
        m_invariant3.reset();
    }

protected:

    std::optional<DeformationGradient> m_deformationGradient;
    std::optional<RightCauchyGreenTensor> m_rightCauchyGreenTensor;

    std::optional<Real> m_determinantF;
    std::optional<GreenLagrangeTensor> m_greenLagrangeTensor;
    std::optional<Real> m_invariant1;
    std::optional<Real> m_invariant2;
    std::optional<Real> m_invariant3;

    void computeDeterminant()
    {
        if (m_deformationGradient.has_value())
        {
            m_determinantF = elasticity::determinantSquareMatrix(*m_deformationGradient);
        }
        else if (m_rightCauchyGreenTensor.has_value())
        {
            m_determinantF = sqrt(elasticity::determinantSquareMatrix(*m_rightCauchyGreenTensor));
        }
        else
        {
            throw std::logic_error("Cannot compute determinant. Neither deformation gradient nor right Cauchy-Green tensor is provided");
        }
    }

    void computeRightCauchyGreenTensor()
    {
        if (m_deformationGradient.has_value())
        {
            const auto& F = *m_deformationGradient;
            m_rightCauchyGreenTensor = F.transposed() * F;
        }
        else
        {
            throw std::logic_error("Cannot compute right Cauchy-Green tensor. Deformation gradient is not provided");
        }
    }

    void computeGreenLagrangeTensor()
    {
        const auto& C = getRightCauchyGreenTensor();
        m_greenLagrangeTensor = static_cast<Real>(0.5) * (C - identity);
    }

    void computeInvariant1()
    {
        const auto& C = getRightCauchyGreenTensor();
        m_invariant1 = sofa::type::trace(C);
    }

    void computeInvariant2()
    {
        const auto I1 = getInvariant1();
        const auto& C = getRightCauchyGreenTensor();
        const auto trC2 = elasticity::squaredFrobeniusNorm(C);
        m_invariant2 = static_cast<Real>(0.5) * (I1 * I1 - trC2);
    }

    void computeInvariant3()
    {
        const auto J = getDeterminantDeformationGradient();
        m_invariant1 = J * J;
    }
};

}

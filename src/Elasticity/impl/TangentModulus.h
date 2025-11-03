#pragma once

#include <sofa/core/trait/DataTypes.h>
#include <sofa/type/MatSym.h>

namespace elasticity
{

template <class DataTypes>
class TangentModulus
{
private:
    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;
    static constexpr sofa::Size spatial_dimension_square = spatial_dimensions * spatial_dimensions;
    using Real = sofa::Real_t<DataTypes>;

public:
    TangentModulus() = delete;

    template<class Callable>
    TangentModulus(Callable callable) : m_matrix(sofa::type::NOINIT)
    {
        fill(callable);
    }

    template<class Callable>
    void fill(Callable callable)
    {
        for (sofa::Size i = 0; i < spatial_dimension_square; ++i)
        {
            const auto a = i / spatial_dimensions;
            const auto b = i % spatial_dimensions;

            for (sofa::Size j = i; j < spatial_dimension_square; ++j) // the reduced representation is symmetric, that is why j starts at i
            {
                const auto c = j / spatial_dimensions;
                const auto d = j % spatial_dimensions;

                m_matrix(i, j) = callable(a, b, c, d);
            }
        }
    }

    Real& operator()(sofa::Size i, sofa::Size j, sofa::Size k, sofa::Size l)
    {
        const auto a = i * spatial_dimensions + j;
        const auto b = k * spatial_dimensions + l;
        return m_matrix(a, b);
    }

    Real operator()(sofa::Size i, sofa::Size j, sofa::Size k, sofa::Size l) const
    {
        const auto a = i * spatial_dimensions + j;
        const auto b = k * spatial_dimensions + l;
        return m_matrix(a, b);
    }

private:
    sofa::type::MatSym<spatial_dimension_square, Real> m_matrix;
};

}

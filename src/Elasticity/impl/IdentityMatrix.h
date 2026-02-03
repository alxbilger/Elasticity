#pragma once

#include <sofa/type/Mat.h>

namespace elasticity
{

struct IdentityMatrix
{};

template<class real>
struct ScaledIdentityMatrix
{
    real scale;
};

template<sofa::Size N, class real>
constexpr sofa::type::Mat<N, N, real> operator+(const elasticity::IdentityMatrix& I, const sofa::type::Mat<N, N, real>& M)
{
    sofa::type::Mat<N, N, real> res(M);
    for (sofa::Size i = 0; i < N; ++i)
    {
        res[i][i] += static_cast<real>(1);
    }
    return res;
}

template<sofa::Size N, class real>
constexpr sofa::type::Mat<N, N, real> operator+(const sofa::type::Mat<N, N, real>& M, const elasticity::IdentityMatrix& I)
{
    return I + M;
}

template<sofa::Size N, class real>
constexpr sofa::type::Mat<N, N, real> operator-(const elasticity::IdentityMatrix& I, const sofa::type::Mat<N, N, real>& M)
{
    sofa::type::Mat<N, N, real> res(-M);
    for (sofa::Size i = 0; i < N; ++i)
    {
        res[i][i] += static_cast<real>(1);
    }
    return res;
}

template<sofa::Size N, class real>
constexpr sofa::type::Mat<N, N, real> operator-(const sofa::type::Mat<N, N, real>& M, const elasticity::IdentityMatrix& I)
{
    sofa::type::Mat<N, N, real> res(M);
    for (sofa::Size i = 0; i < N; ++i)
    {
        res[i][i] -= static_cast<real>(1);
    }
    return res;
}

template<class real>
constexpr ScaledIdentityMatrix<real> operator*(const IdentityMatrix& I, real s)
{
    return { s };
}

template<class real>
constexpr ScaledIdentityMatrix<real> operator*(real s, const IdentityMatrix& I)
{
    return { s };
}

template<sofa::Size N, class real>
constexpr sofa::type::Mat<N, N, real> operator-(const elasticity::ScaledIdentityMatrix<real>& I, const sofa::type::Mat<N, N, real>& M)
{
    sofa::type::Mat<N, N, real> res(-M);
    for (sofa::Size i = 0; i < N; ++i)
    {
        res[i][i] += I.scale;
    }
    return res;
}

template<sofa::Size N, class real>
constexpr sofa::type::Mat<N, N, real> operator-(const sofa::type::Mat<N, N, real>& M, const elasticity::ScaledIdentityMatrix<real>& I)
{
    sofa::type::Mat<N, N, real> res(M);
    for (sofa::Size i = 0; i < N; ++i)
    {
        res[i][i] -= I.scale;
    }
    return res;
}

template<sofa::Size N, class real>
constexpr sofa::type::Mat<N, N, real> operator+(const elasticity::ScaledIdentityMatrix<real>& I, const sofa::type::Mat<N, N, real>& M)
{
    sofa::type::Mat<N, N, real> res(M);
    for (sofa::Size i = 0; i < N; ++i)
    {
        res[i][i] += I.scale;
    }
    return res;
}

template<sofa::Size N, class real>
constexpr sofa::type::Mat<N, N, real> operator+(const sofa::type::Mat<N, N, real>& M, const elasticity::ScaledIdentityMatrix<real>& I)
{
    return I + M;
}



}

#pragma once

#include <sofa/topology/Element.h>
#include <sofa/type/Mat.h>

#include <Eigen/Geometry>
#include <boost/math/policies/policy.hpp>

namespace elasticity
{

template <sofa::Size N, typename real>
real determinantSquareMatrix(const sofa::type::Mat<N, N, real>& mat)
{
    if constexpr (N == 1)
    {
        return mat(0, 0);
    }
    else if constexpr (N == 2)
    {
        return mat(0, 0) * mat(1, 1) - mat(0, 1) * mat(1, 0);
    }
    else
    {
        real det = 0;
        for (size_t p = 0; p < N; ++p)
        {
            sofa::type::Mat<N - 1, N - 1, real> submat;
            for (size_t i = 1; i < N; ++i)
            {
                size_t colIndex = 0;
                for (size_t j = 0; j < N; ++j)
                {
                    if (j == p) continue;
                    submat(i - 1, colIndex++) = mat(i, j);
                }
            }
            det += ((p % 2 == 0) ? 1 : -1) * mat(0, p) * determinantSquareMatrix(submat);
        }
        return det;
    }
}

/**
 * Computes the determinant of a given matrix.
 * For square matrices, the standard determinant is computed.
 * For non-square matrices, the determinant is computed using the product of the transposed
 * matrix and the original matrix, followed by taking the square root of the result.
 *
 * @param mat The input matrix of size LxC.
 * @return The determinant of the matrix.
 */
template <sofa::Size L, sofa::Size C, class real>
real determinant(const sofa::type::Mat<L, C, real>& mat)
{
    if constexpr (L == C)
    {
        return determinantSquareMatrix(mat);
    }
    else
    {
        return std::sqrt(determinantSquareMatrix(mat.transposed() * mat));
    }
}

template <sofa::Size L, sofa::Size C, class real>
sofa::type::Mat<C, L, real> leftPseudoInverse(const sofa::type::Mat<L, C, real>& mat)
{
    return (mat.transposed() * mat).inverted() * mat.transposed();
}

template <sofa::Size L, sofa::Size C, class real>
sofa::type::Mat<C, L, real> rightPseudoInverse(const sofa::type::Mat<L, C, real>& mat)
{
    return mat.transposed() * (mat * mat.transposed()).inverted();
}

/**
 * Computes the inverse of a given matrix.
 * For square matrices (L == C), the standard matrix inverse is computed.
 * For non-square matrices, the left pseudo-inverse is returned.
 *
 * @param mat The input matrix of size LxC to be inverted or pseudo-inverted.
 * @return A matrix of size CxL representing the inverse or left pseudo-inverse of the input matrix.
 */
template <sofa::Size L, sofa::Size C, class real>
sofa::type::Mat<C, L, real> inverse(const sofa::type::Mat<L, C, real>& mat)
{
    if constexpr (L == C)
    {
        return mat.inverted();
    }
    else
    {
        return leftPseudoInverse(mat);
    }
}

/**
 * Stack the columns of a matrix on top of each other to form a vector.
 */
template <sofa::Size L, sofa::Size C, class real>
constexpr auto flatten(const sofa::type::Mat<L, C, real>& mat)
{
    sofa::type::Vec<L * C, real> res(sofa::type::NOINIT);
    sofa::Index index = 0;
    for (sofa::Size c = 0; c < C; ++c)
    {
        for (sofa::Size l = 0; l < L; ++l)
        {
            res[index++] = mat(l, c);
        }
    }
    return res;
}

template <sofa::Size L1, sofa::Size C1, sofa::Size L2, sofa::Size C2, class real>
constexpr sofa::type::Mat<L1 * L2, C1 * C2, real> kroneckerProduct(const sofa::type::Mat<L1, C1, real>& mat1, const sofa::type::Mat<L2, C2, real>& mat2)
{
    sofa::type::Mat<L1 * L2, C1 * C2, real> result(sofa::type::NOINIT);

    for (sofa::Size i1 = 0; i1 < L1; ++i1)
    {
        for (sofa::Size j1 = 0; j1 < C1; ++j1)
        {
            for (sofa::Size i2 = 0; i2 < L2; ++i2)
            {
                for (sofa::Size j2 = 0; j2 < C2; ++j2)
                {
                    result(i1 * L2 + i2, j1 * C2 + j2) = mat1(i1, j1) * mat2(i2, j2);
                }
            }
        }
    }

    return result;
}

template <sofa::Size L1, sofa::Size C1, sofa::Size N, class real>
constexpr sofa::type::Mat<L1 * N, C1, real> kroneckerProduct(const sofa::type::Mat<L1, C1, real>& mat, const sofa::type::Vec<N, real>& vec)
{
    sofa::type::Mat<L1 * N, C1, real> result(sofa::type::NOINIT);

    for (sofa::Size i1 = 0; i1 < L1; ++i1)
    {
        for (sofa::Size j1 = 0; j1 < C1; ++j1)
        {
            for (sofa::Size i2 = 0; i2 < N; ++i2)
            {
                result(i1 * N + i2, j1) = mat(i1, j1) * vec[i2];
            }
        }
    }

    return result;
}

template <sofa::Size N, class real>
constexpr sofa::type::Mat<N, N, real> adjugate(const sofa::type::Mat<N, N, real>& mat)
{
    return elasticity::determinantSquareMatrix(mat) * inverse(mat);
}

template <class real>
constexpr sofa::type::Mat<2, 2, real> adjugate(const sofa::type::Mat<2, 2, real>& mat)
{
    sofa::type::Mat<2, 2, real> adj(sofa::type::NOINIT);
    adj(0, 0) = mat(1, 1);
    adj(0, 1) = -mat(0, 1);
    adj(1, 0) = -mat(1, 0);
    adj(1, 1) = mat(0, 0);
    return adj;
}

template <class real>
constexpr sofa::type::Mat<3, 3, real> adjugate(const sofa::type::Mat<3, 3, real>& mat)
{
    sofa::type::Mat<3, 3, real> adj(sofa::type::NOINIT);
    adj(0, 0) = mat(1, 1) * mat(2, 2) - mat(1, 2) * mat(2, 1);
    adj(0, 1) = mat(0, 2) * mat(2, 1) - mat(0, 1) * mat(2, 2);
    adj(0, 2) = mat(0, 1) * mat(1, 2) - mat(0, 2) * mat(1, 1);
    adj(1, 0) = mat(1, 2) * mat(2, 0) - mat(1, 0) * mat(2, 2);
    adj(1, 1) = mat(0, 0) * mat(2, 2) - mat(0, 2) * mat(2, 0);
    adj(1, 2) = mat(0, 2) * mat(1, 0) - mat(0, 0) * mat(1, 2);
    adj(2, 0) = mat(1, 0) * mat(2, 1) - mat(1, 1) * mat(2, 0);
    adj(2, 1) = mat(0, 1) * mat(2, 0) - mat(0, 0) * mat(2, 1);
    adj(2, 2) = mat(0, 0) * mat(1, 1) - mat(0, 1) * mat(1, 0);
    return adj;
}

template <sofa::Size N, class real>
constexpr real squaredFrobeniusNorm(const sofa::type::Mat<N, N, real>& mat)
{
    real norm = 0;
    for (sofa::Size i = 0; i < N; ++i)
    {
        norm += sofa::type::dot(mat[i], mat[i]);
    }
    return norm;
}

template<sofa::Size N, class real>
struct IdentityMatrix
{};

template<sofa::Size N, class real>
struct ScaledIdentityMatrix
{
    real scale;
};

template<sofa::Size N, class real>
constexpr sofa::type::Mat<N, N, real> operator+(const elasticity::IdentityMatrix<N, real>& I, const sofa::type::Mat<N, N, real>& M)
{
    sofa::type::Mat<N, N, real> res(M);
    for (sofa::Size i = 0; i < N; ++i)
    {
        res[i][i] += static_cast<real>(1);
    }
    return res;
}

template<sofa::Size N, class real>
constexpr sofa::type::Mat<N, N, real> operator+(const sofa::type::Mat<N, N, real>& M, const elasticity::IdentityMatrix<N, real>& I)
{
    sofa::type::Mat<N, N, real> res(M);
    for (sofa::Size i = 0; i < N; ++i)
    {
        res[i][i] += static_cast<real>(1);
    }
    return res;
}

template<sofa::Size N, class real>
constexpr sofa::type::Mat<N, N, real> operator-(const elasticity::IdentityMatrix<N, real>& I, const sofa::type::Mat<N, N, real>& M)
{
    sofa::type::Mat<N, N, real> res(-M);
    for (sofa::Size i = 0; i < N; ++i)
    {
        res[i][i] += static_cast<real>(1);
    }
    return res;
}

template<sofa::Size N, class real>
constexpr sofa::type::Mat<N, N, real> operator-(const sofa::type::Mat<N, N, real>& M, const elasticity::IdentityMatrix<N, real>& I)
{
    sofa::type::Mat<N, N, real> res(M);
    for (sofa::Size i = 0; i < N; ++i)
    {
        res[i][i] -= static_cast<real>(1);
    }
    return res;
}

template<sofa::Size N, class real>
constexpr ScaledIdentityMatrix<N, real> operator*(const IdentityMatrix<N, real>& I, real s)
{
    return { s };
}

template<sofa::Size N, class real>
constexpr ScaledIdentityMatrix<N, real> operator*(real s, const IdentityMatrix<N, real>& I)
{
    return { s };
}

template<sofa::Size N, class real>
constexpr sofa::type::Mat<N, N, real> operator-(const elasticity::ScaledIdentityMatrix<N, real>& I, const sofa::type::Mat<N, N, real>& M)
{
    sofa::type::Mat<N, N, real> res(-M);
    for (sofa::Size i = 0; i < N; ++i)
    {
        res[i][i] += I.scale;
    }
    return res;
}

template<sofa::Size N, class real>
constexpr sofa::type::Mat<N, N, real> operator-(const sofa::type::Mat<N, N, real>& M, const elasticity::ScaledIdentityMatrix<N, real>& I)
{
    sofa::type::Mat<N, N, real> res(M);
    for (sofa::Size i = 0; i < N; ++i)
    {
        res[i][i] -= I.scale;
    }
    return res;
}

template<sofa::Size N, class real>
constexpr sofa::type::Mat<N, N, real> operator+(const elasticity::ScaledIdentityMatrix<N, real>& I, const sofa::type::Mat<N, N, real>& M)
{
    sofa::type::Mat<N, N, real> res(M);
    for (sofa::Size i = 0; i < N; ++i)
    {
        res[i][i] += I.scale;
    }
    return res;
}

template<sofa::Size N, class real>
constexpr sofa::type::Mat<N, N, real> operator+(const sofa::type::Mat<N, N, real>& M, const elasticity::ScaledIdentityMatrix<N, real>& I)
{
    sofa::type::Mat<N, N, real> res(M);
    for (sofa::Size i = 0; i < N; ++i)
    {
        res[i][i] += I.scale;
    }
    return res;
}



}

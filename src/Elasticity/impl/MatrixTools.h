#pragma once

#include <sofa/type/Mat.h>
#include <sofa/topology/Element.h>
#include <Eigen/Geometry>

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

}

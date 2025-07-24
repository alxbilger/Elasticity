#pragma once

#include <sofa/type/Mat.h>
#include <sofa/topology/Element.h>
#include <Eigen/Geometry>

namespace elasticity
{

template <size_t N, typename real>
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

// from
// Müller, Matthias, et al. "A robust method to extract the rotational part of deformations." Proceedings of the 9th International Conference on Motion in Games. 2016.
template<class Real>
void extractRotation(
    const sofa::type::Mat<3,3,Real> &A,
    sofa::type::Quat<Real> &q,
    const unsigned int maxIter)
{

    for (unsigned int iter = 0; iter < maxIter; iter++)
    {
        sofa::type::Mat<3,3,Real> R(sofa::type::NOINIT);
        q.toMatrix(R);
        sofa::type::Vec<3, Real> omega =
            (
                R.col(0).cross(A.col(0)) +
                R.col(1).cross(A.col(1)) +
                R.col(2).cross(A.col(2))
            ) * (1.0 / fabs(
                sofa::type::dot(R.col(0), A.col(0)) +
                sofa::type::dot(R.col(1),
                    A.col(1)) +
                sofa::type::dot(R.col(1), A.col(1))) + 1.0e-9);
        Real w = omega.norm();
        if (w < 1.0e-9)
            break;
        Eigen::Vector3<Real> omegaEigen(omega[0], omega[1], omega[2]);
        Eigen::AngleAxis<Real> angleAxis(w, (1.0 / w) * omegaEigen);
        const auto axis = angleAxis.axis();
        q = sofa::type::Quat<Real>({axis(0), axis(1), axis(2)}, angleAxis.angle()) * q;
        q.normalize();
    }
}


}

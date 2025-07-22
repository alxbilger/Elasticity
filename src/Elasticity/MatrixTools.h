#pragma once

#include <sofa/type/Mat.h>
#include <sofa/topology/Element.h>

namespace elasticity
{

/**
 * Build a matrix where the columns are the coordinates of the nodes in an element
 */
template<class GeometryElement, class VecCoord>
auto nodesMatrix(const sofa::topology::Element<GeometryElement>& nodesIndices, const VecCoord& X)
{
    using Coord = typename VecCoord::value_type;
    static constexpr sofa::Size spatial_dimensions = Coord::static_size;
    static constexpr sofa::Size NumberOfNodesInElement = GeometryElement::NumberOfNodes;
    using Real = typename Coord::value_type;

    sofa::type::Mat<spatial_dimensions, NumberOfNodesInElement, Real> X_element;
    for (sofa::Size i = 0; i < spatial_dimensions; ++i)
    {
        for (sofa::Size j = 0; j < NumberOfNodesInElement; ++j)
        {
            X_element[i][j] = X[nodesIndices[j]][i];
        }
    }

    return X_element;
}

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
sofa::type::Mat<L, C, real> inverse(const sofa::type::Mat<L, C, real>& mat)
{
    if constexpr (L == C)
    {
        return mat.inverted();
    }
    else
    {
        //pseudo inverse
        const sofa::type::Mat<C, C, real> JTJ = mat.transposed() * mat;
        const sofa::type::Mat<C, C, real> JTJ_inverse = JTJ.inverted();
        return mat * JTJ_inverse;
    }
}

}

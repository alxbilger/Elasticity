#pragma once

#include <sofa/type/Vec.h>

namespace elasticity
{

using PairIndices = std::pair<sofa::Size, sofa::Size>;

template<
    sofa::Size Rows,
    sofa::Size Cols,
    sofa::Size NNZ,
    std::array<PairIndices, NNZ> NonZeroEntries,
    class Scalar
>
class StaticSparseMatrix
{
    std::array<Scalar, NNZ> values;

public:
    void set(sofa::Size index, Scalar val)
    {
        values[index] = val;
    }

    Scalar element(sofa::Size i, sofa::Size j) const
    {
        for (sofa::Size k = 0; k < NNZ; ++k)
        {
            auto [ei, ej] = NonZeroEntries[k];
            if (ei == i && ej == j)
            {
                return values[k];
            }
        }
        return 0;
    }

    Scalar& element(sofa::Size i, sofa::Size j)
    {
        for (sofa::Size k = 0; k < NNZ; ++k)
        {
            auto [ei, ej] = NonZeroEntries[k];
            if (ei == i && ej == j)
            {
                return values[k];
            }
        }
        throw std::runtime_error("Element cannot be modified because it has been fixed to zero at compile-time");
    }

    Scalar operator()(sofa::Size i, sofa::Size j) const
    {
        return element(i, j);
    }

    Scalar& operator()(sofa::Size i, sofa::Size j)
    {
        return element(i, j);
    }

    void add(sofa::Size i, sofa::Size j, Scalar val)
    {
        element(i, j) += val;
    }

    void scale(const Scalar scale)
    {
        for (auto& v : values)
            v *= scale;
    }

    sofa::type::Vec<Rows, Scalar> multiply(const sofa::type::Vec<Cols, Scalar>& x) const
    {
        sofa::type::Vec<Rows, Scalar> result{};
        for (sofa::Size k = 0; k < NNZ; ++k)
        {
            auto [i, j] = NonZeroEntries[k];
            result[i] += values[k] * x[j];
        }
        return result;
    }

    sofa::type::Vec<Cols, Scalar> multTranspose(const sofa::type::Vec<Rows, Scalar>& x) const
    {
        sofa::type::Vec<Cols, Scalar> result{};
        for (sofa::Size k = 0; k < NNZ; ++k)
        {
            auto [i, j] = NonZeroEntries[k];
            result[j] += values[k] * x[i];
        }
        return result;
    }

    // template<
    //     sofa::Size OtherMatrixCols,
    //     sofa::Size OtherMatrixNNZ,
    //     std::array<std::pair<sofa::Size, sofa::Size>, NNZ> OtherMatrixNonZeroEntries
    // >
    // auto multiply(const StaticSparseMatrix<Cols, OtherMatrixCols, OtherMatrixNNZ, OtherMatrixNonZeroEntries, Scalar>& matrix)
    // {
    //
    // }
};

template<
    sofa::Size Rows,
    sofa::Size Cols,
    sofa::Size NNZ,
    std::array<PairIndices, NNZ> NonZeroEntries,
    class Scalar
>
sofa::type::Vec<Rows, Scalar> operator*(const StaticSparseMatrix<Rows, Cols, NNZ, NonZeroEntries, Scalar>& matrix, const sofa::type::Vec<Cols, Scalar>& vec)
{
    return matrix.multiply(vec);
}


}

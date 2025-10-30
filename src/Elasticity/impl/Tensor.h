#pragma once
#include <array>

namespace elasticity
{

template<typename T, std::size_t... Dims>
struct Tensor
{
    static constexpr std::size_t rank = sizeof...(Dims);
    static constexpr std::size_t size = (Dims * ...);
    static_assert(size != 0, "None of the dimension of the tensor can be zero.");
    static constexpr std::array dims {Dims...};

private:
    std::array<T, size> data_{};

public:
    // Default constructor
    constexpr Tensor() = default;

    // Access using flat index
    constexpr T& operator[](std::size_t i) { return data_[i]; }
    constexpr const T& operator[](std::size_t i) const { return data_[i]; }

    // Multi-index access
    template<typename... Idx>
    constexpr T& operator()(Idx... idx)
    {
        static_assert(sizeof...(Idx) == rank, "Incorrect number of indices");
        return data_[flattenIndex(idx...)];
    }

    template<typename... Idx>
    constexpr const T& operator()(Idx... idx) const
    {
        static_assert(sizeof...(Idx) == rank, "Incorrect number of indices");
        return data_[flattenIndex(idx...)];
    }

    sofa::type::Mat<dims[0], dims[1], T> toMat() const requires (rank == 2)
    {
        sofa::type::Mat<dims[0], dims[1], T> mat;

        for (std::size_t i = 0; i < dims[0]; ++i)
            for (std::size_t j = 0; j < dims[1]; ++j)
                mat(i, j) = (*this)(i, j);

        return mat;
    }

private:
    // Compute flat index from multi-index at compile-time
    template<typename... Idx>
    static constexpr std::size_t flattenIndex(Idx... idx)
    {
        constexpr std::array<std::size_t, rank> dims = {Dims...};
        std::array<std::size_t, rank> indices = {static_cast<std::size_t>(idx)...};
        std::size_t flat = 0, stride = 1;
        for (std::size_t i = rank; i-- > 0;) {
            flat += indices[i] * stride;
            stride *= dims[i];
        }
        return flat;
    }
};

static_assert(Tensor<double, 1, 1, 1, 1>::rank == 4);
static_assert(Tensor<double, 1, 2, 3, 4>::size == 2*3*4);

}

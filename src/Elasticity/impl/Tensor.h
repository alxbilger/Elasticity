#include <array>

namespace elasticity
{

namespace tensor
{
constexpr struct AllElements{} all;

struct Index
{
    static std::size_t get(AllElements&) { return -1; }
    static std::size_t get(std::size_t i) { return i;}
};
}

template<typename T, std::size_t... Dims>
struct Tensor {
    static constexpr std::size_t rank = sizeof...(Dims);
    static constexpr std::size_t size = (Dims * ...);
    static_assert(size != 0, "None of the dimension of the tensor can be zero.");
    static constexpr std::array dims {Dims...};

    // template<typename TOther, std::size_t... DimsOther>
    // friend struct Tensor<TOther, DimsOther...>;

public:
    std::array<T, size> data_{};

public:
    // Default constructor
    constexpr Tensor() = default;

    // Access using flat index
    constexpr T& operator[](std::size_t i) { return data_[i]; }
    constexpr const T& operator[](std::size_t i) const { return data_[i]; }

    // Multi-index access
    template<typename... Idx>
    constexpr T& operator()(Idx... idx) {
        static_assert(sizeof...(Idx) == rank, "Incorrect number of indices");
        return data_[flatten_index(idx...)];
    }

    template<typename... Idx>
    constexpr const T& operator()(Idx... idx) const {
        static_assert(sizeof...(Idx) == rank, "Incorrect number of indices");
        return data_[flatten_index(idx...)];
    }

    template<class... Indices>
    constexpr auto subtensor(Indices... indices) const
    {
        static_assert(sizeof...(Indices) == rank, "Incorrect number of indices");
        //dimension is not fixed when the corresponding index is not a std::size_t, but a tensor::AllElements
        constexpr std::array isDimensionFixed { !std::is_same_v<tensor::AllElements, std::remove_cvref_t<Indices>> ...};
        static_assert(std::any_of(isDimensionFixed.begin(), isDimensionFixed.end(), [](auto x) { return x; }), "At least one dimension must be fixed");

        const std::array<std::size_t, rank> integerIndices { tensor::Index::get(indices)...};

        constexpr auto nbFixedDimensions = std::count(isDimensionFixed.begin(), isDimensionFixed.end(), true);
        constexpr auto rankSubTensor = rank - nbFixedDimensions;

        constexpr auto dimensionSubTensor = []()
        {
            std::array<std::size_t, rankSubTensor> out{};
            std::size_t i = 0;
            for (std::size_t j = 0; j < rank; ++j)
            {
                if (isDimensionFixed[j])
                    continue;
                out[i++] = dims[j];
            }
            return out;
        }();

        auto subTensor = [&]<std::size_t... I>(std::index_sequence<I...>) -> Tensor<T, dimensionSubTensor[I]...> {
            return {};
        }(std::make_index_sequence<rankSubTensor>{});

        for (std::size_t i = 0, j = 0; i < size; ++i)
        {
            auto ids = unflatten(i);
            bool isIn = true;
            for (std::size_t k = 0; k < rank; ++k)
            {
                if (isDimensionFixed[k])
                {
                    if (integerIndices[k] != ids[k])
                    {
                        isIn = false;
                        break;
                    }
                }
            }
            if (isIn)
            {
                subTensor.data_[j++] = data_[i];
            }
        }

        return subTensor;
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
    static constexpr std::size_t flatten_index(Idx... idx) {
        constexpr std::array<std::size_t, rank> dims = {Dims...};
        std::array<std::size_t, rank> indices = {static_cast<std::size_t>(idx)...};
        std::size_t flat = 0, stride = 1;
        for (std::size_t i = rank; i-- > 0;) {
            flat += indices[i] * stride;
            stride *= dims[i];
        }
        return flat;
    }

    static constexpr std::array<std::size_t, rank> unflatten(std::size_t i)
    {
        std::array<std::size_t, rank> idx{};
        for (std::size_t d = rank; d-- > 0;) {
            idx[d] = i % dims[d];
            i /= dims[d];
        }
        return idx;
    }
};

static_assert(Tensor<double, 1, 1, 1, 1>::rank == 4);
static_assert(Tensor<double, 1, 2, 3, 4>::size == 2*3*4);
static_assert(std::same_as<decltype(Tensor<double, 6, 2, 3, 32>{}.subtensor(tensor::all, 2, tensor::all, 3)), Tensor<double, 6, 3>>);
static_assert(std::same_as<decltype(Tensor<double, 6, 2, 3, 32>{}.subtensor(4, 2, tensor::all, tensor::all)), Tensor<double, 3, 32>>);

}

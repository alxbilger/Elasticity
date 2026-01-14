#pragma once

#include <Elasticity/impl/VecSoA.h>
#include <sofa/type/Mat.h>

namespace elasticity
{

template <sofa::Size L, sofa::Size C, class real>
class MatSoA
{
private:
    std::array<sofa::type::vector<real>, L * C> m_data;

public:
    MatSoA() = default;
    explicit MatSoA(std::size_t count) { resize(count); }

    void resize(std::size_t size)
    {
        if (this->size() != size)
        {
            std::for_each(m_data.begin(), m_data.end(), [size](std::vector<real>& vec){ vec.resize(size); });
        }
    }

    std::size_t size() const { return m_data.front().size(); }

    sofa::type::vector<real>& element(std::size_t row, std::size_t col)
    {
        return m_data[row * C + col];
    }

    const sofa::type::vector<real>& element(std::size_t row, std::size_t col) const
    {
        return m_data[row * C + col];
    }

    void toAoS(sofa::type::vector<sofa::type::Mat<L, C, real> >& aos)
    {
        const auto count = this->size();
        aos.resize(count);

        for (sofa::Size i = 0; i < L; ++i)
        {
            for (sofa::Size j = 0; j < C; ++j)
            {
                for (std::size_t k = 0; k < count; ++k)
                {
                    const auto& data = this->element(i, j);
                    aos[k][i][j] = data[k];
                }
            }
        }
    }
};

template<sofa::Size L, sofa::Size C, class real>
void matrixVectorProduct(VecSoA<L, real>& result, const MatSoA<L, C, real>& mat, const VecSoA<C, real>& vec)
{
    const size_t n = vec.size();
    result.resize(n);

    for (std::size_t row = 0; row < L; ++row)
    {
        auto& result_comp = result.element(row);

        std::fill(result_comp.begin(), result_comp.end(), real(0));

        for (std::size_t col = 0; col < C; ++col)
        {
            const auto& mat_elem = mat.element(row, col);
            const auto& vec_comp = vec.element(col);

            for (size_t i = 0; i < n; ++i)
            {
                result_comp[i] += mat_elem[i] * vec_comp[i];
            }
        }
    }
}



}

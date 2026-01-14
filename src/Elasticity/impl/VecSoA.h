#pragma once

#include <sofa/type/vector.h>
#include <sofa/type/Vec.h>
#include <algorithm>
#include <functional>

namespace elasticity
{

template <sofa::Size N, typename real>
class VecSoA
{
private:
    std::array<sofa::type::vector<real>, N> m_data;

public:

    VecSoA() = default;
    explicit VecSoA(std::size_t count) { resize(count); }

    void resize(std::size_t size)
    {
        if (this->size() != size)
        {
            std::for_each(m_data.begin(), m_data.end(), [size](std::vector<real>& vec){ vec.resize(size); });
        }
    }

    std::size_t size() const { return m_data.front().size(); }

    sofa::type::vector<real>& element(std::size_t i)
    {
        return m_data[i];
    }

    const sofa::type::vector<real>& element(std::size_t i) const
    {
        return m_data[i];
    }

    template<class BinaryOp>
    static void BinaryOperation(VecSoA& result, const VecSoA& a, const VecSoA& b, const BinaryOp& op)
    {
        assert(a.size() == b.size());
        const auto commonSize = a.size();
        result.resize(commonSize);

        for (sofa::Size i = 0; i < N; ++i)
        {
            const auto& a_i = a.m_data[i];
            const auto& b_i = b.m_data[i];
            auto& result_i = result.m_data[i];
            for (std::size_t j = 0; j < commonSize; ++j)
            {
                result_i[j] = op(a_i[j], b_i[j]);
            }
        }
    }

    static void Add(VecSoA& result, const VecSoA& a, const VecSoA& b)
    {
        assert(a.size() == b.size());
        const auto commonSize = a.size();
        result.resize(commonSize);

        for (sofa::Size i = 0; i < N; ++i)
        {
            const auto& a_i = a.m_data[i];
            const auto& b_i = b.m_data[i];
            auto& result_i = result.m_data[i];
            for (std::size_t j = 0; j < commonSize; ++j)
            {
                result_i[j] = a_i[j] + b_i[j];
            }
        }
    }

    void toAoS(sofa::type::vector<sofa::type::Vec<N, real> >& aos)
    {
        const auto count = this->size();
        aos.resize(count);

        for (sofa::Size i = 0; i < N; ++i)
        {
            for (std::size_t j = 0; j < count; ++j)
            {
                aos[j][i] = m_data[i][j];
            }
        }
    }
};

}

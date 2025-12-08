#pragma once
#include <sofa/core/trait/DataTypes.h>
#include <sofa/type/Mat.h>

namespace elasticity
{

template<class DataType>
struct IsotropicElasticityTensor
{
    static constexpr auto spatial_dimensions = DataType::spatial_dimensions;
    static constexpr auto NbIndependentElements = symmetric_tensor::NumberOfIndependentElements<spatial_dimensions>;

    explicit IsotropicElasticityTensor(const sofa::type::Mat<NbIndependentElements, NbIndependentElements, sofa::Real_t<DataType>>& mat) : C(mat) {}
    IsotropicElasticityTensor() = default;



    sofa::type::Vec<NbIndependentElements, sofa::Real_t<DataType>> operator*(const sofa::type::Vec<NbIndependentElements, sofa::Real_t<DataType>>& v) const
    {
        return C * v;
    }

    template<sofa::Size C>
    sofa::type::Mat<NbIndependentElements, C, sofa::Real_t<DataType>> operator*(const sofa::type::Mat<NbIndependentElements, C, sofa::Real_t<DataType>>& v) const noexcept
    {
        return C * v;
    }

    sofa::Real_t<DataType> operator()(sofa::Size i, sofa::Size j) const
    {
        return C(i, j);
    }

    sofa::Real_t<DataType>& operator()(sofa::Size i, sofa::Size j)
    {
        return C(i, j);
    }

    const sofa::type::Mat<NbIndependentElements, NbIndependentElements, sofa::Real_t<DataType>>& toMat() const
    {
        return C;
    }

private:
    sofa::type::Mat<NbIndependentElements, NbIndependentElements, sofa::Real_t<DataType>> C;

};

}

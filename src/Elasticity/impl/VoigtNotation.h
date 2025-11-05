#pragma once

#include <sofa/defaulttype/VecTypes.h>

#include <array>

namespace elasticity
{

template<class DataTypes>
constexpr auto voigtIndices(std::size_t i)
{
    assert(i < symmetric_tensor::NumberOfIndependentElements<DataTypes::spatial_dimensions>);
    if constexpr (DataTypes::spatial_dimensions == 3)
    {
        static constexpr std::array voigt3d {
            std::make_pair(0, 0),
            std::make_pair(1, 1),
            std::make_pair(2, 2),
            std::make_pair(1, 2),
            std::make_pair(0, 2),
            std::make_pair(0, 1)
        };
        assert(i < voigt3d.size());
        return voigt3d[i];
    }
    else if constexpr (DataTypes::spatial_dimensions == 2)
    {
        static constexpr std::array voigt2d {
            std::make_pair(0, 0),
            std::make_pair(1, 1),
            std::make_pair(0, 1)
        };
        assert(i < voigt2d.size());
        return voigt2d[i];
    }
    else
    {
        return std::make_pair(0, 0);
    }
}

template<class DataTypes>
constexpr std::size_t voigtIndex(std::size_t i, std::size_t j)
{
    assert(i < DataTypes::spatial_dimensions);
    assert(j < DataTypes::spatial_dimensions);
    if (i == j)
        return i;
    return symmetric_tensor::NumberOfIndependentElements<DataTypes::spatial_dimensions> - i - j;
}

static_assert(voigtIndex<sofa::defaulttype::Vec3Types>(0,0) == 0);
static_assert(voigtIndex<sofa::defaulttype::Vec3Types>(0,1) == 5);
static_assert(voigtIndex<sofa::defaulttype::Vec3Types>(0,2) == 4);
static_assert(voigtIndex<sofa::defaulttype::Vec3Types>(1,0) == 5);
static_assert(voigtIndex<sofa::defaulttype::Vec3Types>(1,1) == 1);
static_assert(voigtIndex<sofa::defaulttype::Vec3Types>(1,2) == 3);
static_assert(voigtIndex<sofa::defaulttype::Vec3Types>(2,0) == 4);
static_assert(voigtIndex<sofa::defaulttype::Vec3Types>(2,1) == 3);
static_assert(voigtIndex<sofa::defaulttype::Vec3Types>(2,2) == 2);

static_assert(voigtIndex<sofa::defaulttype::Vec2Types>(0,0) == 0);
static_assert(voigtIndex<sofa::defaulttype::Vec2Types>(0,1) == 2);
static_assert(voigtIndex<sofa::defaulttype::Vec2Types>(1,0) == 2);
static_assert(voigtIndex<sofa::defaulttype::Vec2Types>(1,1) == 1);

static_assert(voigtIndex<sofa::defaulttype::Vec1Types>(0,0) == 0);


}

#pragma once

#include <sofa/core/objectmodel/BaseObject.h>
#include <variant>

namespace elasticity
{


template <class DataTypes, class ElementType, class... Methods>
struct RotationMethodsContainer
{
    std::variant<Methods...> m_rotationComputer;

    explicit RotationMethodsContainer(sofa::core::objectmodel::BaseObject* parent)
        : d_rotationMethod(parent->initData(&d_rotationMethod, "rotationMethod", ("The method used to compute the element rotations.\n" + RotationMethodsItems::dataDescription()).c_str()))
    {}

    using RotationMatrix = sofa::type::Mat<DataTypes::spatial_dimensions, DataTypes::spatial_dimensions, sofa::Real_t<DataTypes>>;
    static constexpr sofa::Size NumberOfNodesInElement = ElementType::NumberOfNodes;

    void computeRotation(RotationMatrix& rotationMatrix, RotationMatrix& initialRotationMatrix,
        const std::array<sofa::Coord_t<DataTypes>, NumberOfNodesInElement>& nodesPosition,
        const std::array<sofa::Coord_t<DataTypes>, NumberOfNodesInElement>& nodesRestPosition)
    {
        std::visit(
            [&](auto& rotationComputer)
            {
                rotationComputer.computeRotation(rotationMatrix, initialRotationMatrix, nodesPosition, nodesRestPosition);
            },
            m_rotationComputer);
    }

    MAKE_SELECTABLE_ITEMS(RotationMethodsItems, Methods::getItem()...);
    sofa::Data< RotationMethodsItems > d_rotationMethod;

    void selectRotationMethod()
    {
        const std::size_t selectedId = d_rotationMethod.getValue();

        if (selectedId < std::variant_size_v<decltype(m_rotationComputer)>)
        {
            doSelectRotationMethod<std::variant_size_v<decltype(m_rotationComputer)> - 1>(selectedId);
        }
    }

private:

    template<std::size_t Id>
    void doSelectRotationMethod(const std::size_t selectedId)
    {
        if (selectedId == Id)
        {
            m_rotationComputer.template emplace<std::variant_alternative_t<Id, decltype(m_rotationComputer)>>();
        }
        else
        {
            if constexpr (Id > 0)
            {
                doSelectRotationMethod<Id - 1>(selectedId);
            }
        }
    }
};

}

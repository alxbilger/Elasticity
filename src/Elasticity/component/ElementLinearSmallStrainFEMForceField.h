#pragma once
#include <Elasticity/component/LinearSmallStrainFEMForceField.h>

namespace elasticity
{

template <class DataTypes, class ElementType>
class ElementLinearSmallStrainFEMForceField : public LinearSmallStrainFEMForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE2(ElementLinearSmallStrainFEMForceField, DataTypes, ElementType),
               SOFA_TEMPLATE(LinearSmallStrainFEMForceField, DataTypes));

    static const std::string GetCustomClassName()
    {
        return std::string(sofa::geometry::elementTypeToString(ElementType::Element_type)) +
               "LinearSmallStrainFEMForceField";
    }

    static const std::string GetCustomTemplateName() { return DataTypes::Name(); }

protected:
    void selectFEMTypes() override;
};

}  // namespace elasticity

#pragma once
#include <Elasticity/component/LinearSmallStrainFEMForceField.h>

namespace elasticity
{

/**
 * This class is a @LinearSmallStrainFEMForceField, but working on a single specific type of
 * element.
 *
 * While @LinearSmallStrainFEMForceField works on all the supported in the topology, this class
 * works only a given element type, specified at compile-time.
 *
 * The component name is customized to include the element type in its name.
 */
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

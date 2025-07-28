#pragma once

#include <Elasticity/component/LinearSmallStrainFEMForceField.h>

namespace elasticity
{
template <class DataTypes>
class CorotationalFEMForceField : public LinearSmallStrainFEMForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(CorotationalFEMForceField, DataTypes),
               SOFA_TEMPLATE(sofa::core::behavior::ForceField, DataTypes));

protected:
    void selectFEMTypes() override;

private:
    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;

    template<class ElementType>
    void addCorotationalFEMType();
};

}  // namespace elasticity

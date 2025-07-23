#pragma once

#include <Elasticity/LinearSmallStrainFEMForceField.h>

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
};

}  // namespace elasticity

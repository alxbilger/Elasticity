#pragma once
#include <Elasticity/config.h>
#include <sofa/core/behavior/ForceField.h>

namespace elasticity
{
template <class DataTypes, class ElementType>
class SoAElementLinearSmallStrainFEMForceField : public sofa::core::behavior::ForceField<DataTypes>
{
public:
    SOFA_CLASS(SoAElementLinearSmallStrainFEMForceField,
           sofa::core::behavior::ForceField<DataTypes>);

    static constexpr auto spatial_dimensions = DataTypes::spatial_dimensions;

    void addForce(const sofa::core::MechanicalParams* mparams,
                  sofa::DataVecDeriv_t<DataTypes>& f,
                  const sofa::DataVecCoord_t<DataTypes>& x,
                  const sofa::DataVecDeriv_t<DataTypes>& v) override;

    void addDForce(const sofa::core::MechanicalParams* mparams,
                   sofa::DataVecDeriv_t<DataTypes>& df,
                   const sofa::DataVecDeriv_t<DataTypes>& dx) override;

    SReal getPotentialEnergy(
        const sofa::core::MechanicalParams*,
        const sofa::DataVecDeriv_t<DataTypes>& x) const override;

};

}  // namespace elasticity

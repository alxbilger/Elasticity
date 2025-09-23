#pragma once

#include <Elasticity/component/BaseElasticityFEMForceField.h>

namespace elasticity
{

template <class DataTypes>
class HyperelasticityFEMForceField : public BaseElasticityFEMForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(HyperelasticityFEMForceField, DataTypes),
               SOFA_TEMPLATE(BaseElasticityFEMForceField, DataTypes));

private:
    using DataVecCoord = sofa::DataVecDeriv_t<DataTypes>;
    using DataVecDeriv = sofa::DataVecDeriv_t<DataTypes>;
    using Real = sofa::Real_t<DataTypes>;

public:

    void addForce(const sofa::core::MechanicalParams* mparams, DataVecDeriv& f,
                  const DataVecCoord& x, const DataVecDeriv& v) override;

    void addDForce(const sofa::core::MechanicalParams* mparams, DataVecDeriv& df,
                   const DataVecDeriv& dx) override;

    void buildStiffnessMatrix(sofa::core::behavior::StiffnessMatrix* matrix) override;

    SReal getPotentialEnergy(const sofa::core::MechanicalParams*,
                             const DataVecCoord& x) const override;

protected:
    void selectFEMTypes() override;
};

}  // namespace elasticity

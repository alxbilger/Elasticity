#pragma once

#include <Elasticity/component/BaseElasticityFEMForceField.h>

#include <Elasticity/impl/FEM/NonLinearFEM.h>

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

    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;

public:

    SReal getPotentialEnergy(const sofa::core::MechanicalParams*,
                             const DataVecCoord& x) const override;

protected:
    void selectFEMTypes() override;

    template<class ElementType>
    void addNonLinearFEMType();

    void applyLambda(const std::function<void(BaseFEM<DataTypes>&)>& callable) override;

    sofa::type::vector<std::unique_ptr<BaseNonLinearFEM<DataTypes>>> m_finiteElements;
};

}  // namespace elasticity

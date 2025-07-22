#pragma once

#include <Elasticity/BaseLinearSmallStrainFEMForceField.h>
#include <sofa/core/topology/BaseMeshTopology.h>

#include <Elasticity/LinearFEM.h>

namespace elasticity
{

template <class DataTypes>
class LinearSmallStrainFEMForceField : public BaseLinearSmallStrainFEMForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(LinearSmallStrainFEMForceField, DataTypes),
               SOFA_TEMPLATE(BaseLinearSmallStrainFEMForceField, DataTypes));

private:
    using DataVecCoord = sofa::DataVecDeriv_t<DataTypes>;
    using DataVecDeriv = sofa::DataVecDeriv_t<DataTypes>;
    using Real = sofa::Real_t<DataTypes>;

    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;

public:
    using BaseLinearSmallStrainFEMForceField<DataTypes>::l_topology;
    using BaseLinearSmallStrainFEMForceField<DataTypes>::d_poissonRatio;
    using BaseLinearSmallStrainFEMForceField<DataTypes>::d_youngModulus;

    void init() override;

    void addForce(const sofa::core::MechanicalParams* mparams, DataVecDeriv& f,
                  const DataVecCoord& x, const DataVecDeriv& v) override;

    void addDForce(const sofa::core::MechanicalParams* mparams, DataVecDeriv& df,
                   const DataVecDeriv& dx) override;

    void buildStiffnessMatrix(sofa::core::behavior::StiffnessMatrix* matrix) override;

    SReal getPotentialEnergy(const sofa::core::MechanicalParams*,
                             const DataVecCoord& x) const override;

protected:
    void precomputeElementStiffness() override;

    virtual void selectFEMTypes();

    template<class ElementType>
    void addFEMType();

    sofa::type::vector<std::unique_ptr<BaseLinearFEM<DataTypes>>> m_finiteElements;
};

}  // namespace elasticity

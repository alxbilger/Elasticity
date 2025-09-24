#pragma once

#include <Elasticity/impl/FEM/LinearFEM.h>
#include <Elasticity/component/BaseElasticityFEMForceField.h>
#include <sofa/core/topology/BaseMeshTopology.h>

#include <Elasticity/impl/VonMisesStressContainer.h>

namespace elasticity
{

template <class DataTypes>
class LinearSmallStrainFEMForceField : public BaseElasticityFEMForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(LinearSmallStrainFEMForceField, DataTypes),
               SOFA_TEMPLATE(BaseElasticityFEMForceField, DataTypes));

private:
    using DataVecCoord = sofa::DataVecDeriv_t<DataTypes>;
    using DataVecDeriv = sofa::DataVecDeriv_t<DataTypes>;
    using Real = sofa::Real_t<DataTypes>;

    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;

public:

    sofa::Data<Real> d_poissonRatio;
    sofa::Data<Real> d_youngModulus;
    sofa::Data<bool> d_computeVonMisesStress;
    sofa::Data<sofa::type::vector<Real> > d_vonMisesStressValues;

    void init() override;


    void addForce(const sofa::core::MechanicalParams* mparams, DataVecDeriv& f,
                  const DataVecCoord& x, const DataVecDeriv& v) override;

    SReal getPotentialEnergy(const sofa::core::MechanicalParams*,
                             const DataVecCoord& x) const override;

protected:

    LinearSmallStrainFEMForceField();

    /**
     * With linear small strain, the element stiffness matrix is constant, so it can be precomputed.
     */
    void precomputeElementStiffness();

    void selectFEMTypes() override;

    template<class ElementType>
    void addLinearFEMType();

    sofa::type::vector<std::unique_ptr<BaseLinearFEM<DataTypes>>> m_finiteElements;

    VonMisesStressContainer<Real> m_vonMisesStressContainer;

    void computeVonMisesStress(const DataVecCoord& x);

    void applyLambda(const std::function<void(BaseFEM<DataTypes>&)>& callable) override;
};

}  // namespace elasticity

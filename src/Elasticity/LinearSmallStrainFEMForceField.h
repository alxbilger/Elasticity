#pragma once

#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/topology/BaseMeshTopology.h>

#include <Elasticity/LinearFEM.h>

namespace elasticity
{

template <class DataTypes>
class LinearSmallStrainFEMForceField : public sofa::core::behavior::ForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(LinearSmallStrainFEMForceField, DataTypes),
               SOFA_TEMPLATE(sofa::core::behavior::ForceField, DataTypes));

private:
    using DataVecCoord = sofa::DataVecDeriv_t<DataTypes>;
    using DataVecDeriv = sofa::DataVecDeriv_t<DataTypes>;
    using Real = sofa::Real_t<DataTypes>;

    static constexpr sofa::Size spatial_dimensions = DataTypes::spatial_dimensions;

public:
    /// The topology will give access to the elements
    sofa::SingleLink<LinearSmallStrainFEMForceField, sofa::core::topology::BaseMeshTopology,
                     sofa::BaseLink::FLAG_STOREPATH | sofa::BaseLink::FLAG_STRONGLINK> l_topology;

    sofa::Data<Real> d_poissonRatio;
    sofa::Data<Real> d_youngModulus;

    void init() override;

    void addForce(const sofa::core::MechanicalParams* mparams, DataVecDeriv& f,
                  const DataVecCoord& x, const DataVecDeriv& v) override;

    void addDForce(const sofa::core::MechanicalParams* mparams, DataVecDeriv& df,
                   const DataVecDeriv& dx) override;

    void buildStiffnessMatrix(sofa::core::behavior::StiffnessMatrix* matrix) override;

    SReal getPotentialEnergy(const sofa::core::MechanicalParams*,
                             const DataVecCoord& x) const override;

protected:

    LinearSmallStrainFEMForceField();

    /**
     * Ensure a link to a valid topology. Without a topology, this force field cannot have access
     * to the list of elements.
     */
    void validateTopology();

    /**
     * With linear small strain, the element stiffness matrix is constant, so it can be precomputed.
     */
    void precomputeElementStiffness();

    virtual void selectFEMTypes();

    template<class ElementType>
    void addFEMType();

    sofa::type::vector<std::unique_ptr<BaseLinearFEM<DataTypes>>> m_finiteElements;
};

}  // namespace elasticity

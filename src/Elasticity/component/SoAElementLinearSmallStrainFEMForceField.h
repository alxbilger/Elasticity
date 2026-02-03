#pragma once
#include <Elasticity/config.h>
#include <Elasticity/impl/MatSoA.h>
#include <Elasticity/impl/VecSoA.h>
#include <Elasticity/impl/trait.h>
#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/behavior/TopologyAccessor.h>

namespace elasticity
{
template <class DataTypes, class ElementType>
class SoAElementLinearSmallStrainFEMForceField :
    public sofa::core::behavior::ForceField<DataTypes>,
    public sofa::core::behavior::TopologyAccessor
{
public:
    SOFA_CLASS2(SoAElementLinearSmallStrainFEMForceField,
           sofa::core::behavior::ForceField<DataTypes>, sofa::core::behavior::TopologyAccessor);

    static constexpr auto spatial_dimensions = DataTypes::spatial_dimensions;

    void init() override;

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

private:

    using trait = elasticity::trait<DataTypes, ElementType>;

    VecSoA<trait::NumberOfDofsInElement, sofa::Real_t<DataTypes> > m_elementDisplacement;
    MatSoA<trait::NumberOfDofsInElement, trait::NumberOfDofsInElement, sofa::Real_t<DataTypes> > m_elementStiffness;

    VecSoA<trait::NumberOfDofsInElement, sofa::Real_t<DataTypes> > m_elementForce;

};

}  // namespace elasticity

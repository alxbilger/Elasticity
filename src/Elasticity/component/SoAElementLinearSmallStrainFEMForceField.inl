#include <Elasticity/component/SoAElementLinearSmallStrainFEMForceField.h>
#include <Elasticity/impl/MatSoA.h>
#include <Elasticity/impl/VecSoA.h>

namespace elasticity
{

template <class DataTypes, class ElementType>
void SoAElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::addForce(
    const sofa::core::MechanicalParams*,
    sofa::DataVecDeriv_t<DataTypes>& f,
    const sofa::DataVecCoord_t<DataTypes>& x,
    const sofa::DataVecDeriv_t<DataTypes>& v)
{
    VecSoA<spatial_dimensions, sofa::Real_t<DataTypes> > displacement;
    MatSoA<spatial_dimensions, spatial_dimensions, sofa::Real_t<DataTypes> > stiffness;

    VecSoA<spatial_dimensions, sofa::Real_t<DataTypes> > elementForce;

    elementForce = stiffness * displacement;
}

template <class DataTypes, class ElementType>
void SoAElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::addDForce(
    const sofa::core::MechanicalParams* mparams, sofa::DataVecDeriv_t<DataTypes>& df,
    const sofa::DataVecDeriv_t<DataTypes>& dx)
{
}

template <class DataTypes, class ElementType>
SReal SoAElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::getPotentialEnergy(
    const sofa::core::MechanicalParams*,
    const sofa::DataVecDeriv_t<DataTypes>& x) const
{
    return 0;
}

}  // namespace elasticity

#pragma once

#include <Elasticity/component/material/OgdenMaterial.h>
#include <Elasticity/component/PK2HyperelasticMaterial.inl>

namespace elasticity
{

template <class DataTypes>
OgdenMaterial<DataTypes>::OgdenMaterial()
    : m_mu(initData(&m_mu, static_cast<Real>(1e3), "mu",
                      "Material constant relevant to the shear modulus"))
    , m_alpha(initData(&m_alpha, static_cast<Real>(1e3), "alpha",
                      "Material constant exponent related to strain-stiffening"))
    , m_kappa(initData(&m_kappa, static_cast<Real>(1e3), "kappa",
                      "Material constant related to the bulk modulus"))
{}

template <class DataTypes>
auto OgdenMaterial<DataTypes>::secondPiolaKirchhoffStress(Strain<DataTypes>& strain) -> StressTensor
{
    return sofa::type::Mat<spatial_dimensions,spatial_dimensions,Real>();
}

template <class DataTypes>
auto OgdenMaterial<DataTypes>::elasticityTensor(Strain<DataTypes>& strain) -> ElasticityTensor
{
    return elasticity::FullySymmetric4Tensor<sofa::defaulttype::StdVectorTypes<DataTypes::Coord,DataTypes::Coord,Real>>();
}

}  // namespace elasticity

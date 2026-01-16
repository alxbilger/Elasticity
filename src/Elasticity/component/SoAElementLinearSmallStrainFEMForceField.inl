#pragma once
#include <Elasticity/component/SoAElementLinearSmallStrainFEMForceField.h>


namespace elasticity
{

template <class DataTypes, class ElementType>
void SoAElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::init()
{
    Inherit1::init();
    Inherit2::init();

    m_elementStiffness.resize(l_topology->getNbElements<ElementType>());
    m_elementDisplacement.resize(l_topology->getNbElements<ElementType>());
    m_elementForce.resize(l_topology->getNbElements<ElementType>());
}

template <class DataTypes, class ElementType>
void SoAElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::addForce(
    const sofa::core::MechanicalParams*,
    sofa::DataVecDeriv_t<DataTypes>& f,
    const sofa::DataVecCoord_t<DataTypes>& x,
    const sofa::DataVecDeriv_t<DataTypes>& v)
{
    const auto elements = trait::FiniteElement::getElementSequence(*this->l_topology);
    auto positionAccessor = sofa::helper::getReadAccessor(x);
    auto forceAccessor = sofa::helper::getWriteOnlyAccessor(f);
    auto restPositionAccessor = this->mstate->readRestPositions();

    for (sofa::Size j = 0; j < trait::NumberOfNodesInElement; ++j)
    {
        for (sofa::Size dim = 0; dim < trait::spatial_dimensions; ++dim)
        {
            for (std::size_t elementId = 0; elementId < elements.size(); ++elementId)
            {
                const auto& element = elements[elementId];
                const auto nodeId = element[j];

                m_elementDisplacement.element(j * trait::spatial_dimensions + dim)[elementId] =
                    positionAccessor[nodeId][dim] - restPositionAccessor[nodeId][dim];
            }
        }
    }

    m_elementForce = m_elementStiffness * m_elementDisplacement;

    for (std::size_t elementId = 0; elementId < elements.size(); ++elementId)
    {
        const auto& element = elements[elementId];
        for (sofa::Size j = 0; j < trait::NumberOfNodesInElement; ++j)
        {
            auto& nodeForce = forceAccessor[element[j]];
            for (sofa::Size k = 0; k < trait::spatial_dimensions; ++k)
            {
                // nodeForce[k] -= m_elementForce[j * trait::spatial_dimensions + k];
            }
        }
    }
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

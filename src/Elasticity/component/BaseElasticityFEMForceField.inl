#pragma once

#include <elasticity/component/BaseElasticityFEMForceField.h>
#include <sofa/core/ObjectFactory.h>

namespace elasticity
{

template <class DataTypes>
BaseElasticityFEMForceField<DataTypes>::BaseElasticityFEMForceField()
    : l_topology(initLink("topology", "Link to a topology containing elements"))
{
}

template <class DataTypes>
void BaseElasticityFEMForceField<DataTypes>::init()
{
    sofa::core::behavior::ForceField<DataTypes>::init();

    if (!this->isComponentStateInvalid())
    {
        this->validateTopology();
    }

    if (!this->isComponentStateInvalid())
    {
        selectFEMTypes();
    }
}

template <class DataTypes>
void BaseElasticityFEMForceField<DataTypes>::addForce(const sofa::core::MechanicalParams* mparams,
                                                      DataVecDeriv& f, const DataVecCoord& x,
                                                      const DataVecDeriv& v)
{
    SOFA_UNUSED(mparams);
    SOFA_UNUSED(v);

    auto forceAccessor = sofa::helper::getWriteOnlyAccessor(f);
    auto positionAccessor = sofa::helper::getReadAccessor(x);
    auto restPositionAccessor = this->mstate->readRestPositions();

    applyLambda([&forceAccessor, &positionAccessor, &restPositionAccessor]
        (BaseFEM<DataTypes>& finiteElement)
        {
            finiteElement.addForce(forceAccessor.wref(), positionAccessor.ref(), restPositionAccessor.ref());
        });
}

template <class DataTypes>
void BaseElasticityFEMForceField<DataTypes>::addDForce(const sofa::core::MechanicalParams* mparams,
                                                       DataVecDeriv& df, const DataVecDeriv& dx)
{
    auto dfAccessor = sofa::helper::getWriteAccessor(df);
    auto dxAccessor = sofa::helper::getReadAccessor(dx);
    dfAccessor.resize(dxAccessor.size());

    const Real kFactor = (Real)sofa::core::mechanicalparams::kFactorIncludingRayleighDamping(
        mparams, this->rayleighStiffness.getValue());

    applyLambda([&dfAccessor, &dxAccessor, kFactor]
        (BaseFEM<DataTypes>& finiteElement)
        {
            finiteElement.addDForce(dfAccessor.wref(), dxAccessor.ref(), kFactor);
        });
}
template <class DataTypes>
void BaseElasticityFEMForceField<DataTypes>::buildStiffnessMatrix(
    sofa::core::behavior::StiffnessMatrix* matrix)
{
    auto dfdx = matrix->getForceDerivativeIn(this->mstate)
       .withRespectToPositionsIn(this->mstate);

    applyLambda([&dfdx]
        (BaseFEM<DataTypes>& finiteElement)
        {
            finiteElement.buildStiffnessMatrix(dfdx);
        });
}

template <class DataTypes>
void BaseElasticityFEMForceField<DataTypes>::validateTopology()
{
    if (l_topology.empty())
    {
        msg_info() << "Link to Topology container should be set to ensure right behavior. First "
                      "Topology found in current context will be used.";
        l_topology.set(this->getContext()->getMeshTopologyLink());
    }

    if (l_topology == nullptr)
    {
        msg_error() << "No topology component found at path: " << this->l_topology.getLinkedPath()
                    << ", nor in current context: " << this->getContext()->name
                    << ". Object must have a BaseMeshTopology. "
                    << "The list of available BaseMeshTopology components is: "
                    << sofa::core::ObjectFactory::getInstance()
                           ->listClassesDerivedFrom<sofa::core::topology::BaseMeshTopology>();
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
    }
}

}  // namespace elasticity

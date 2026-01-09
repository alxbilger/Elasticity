#pragma once
#include <Elasticity/impl/ComputeStrategy.h>
#include <sofa/core/visual/DrawMesh.h>
#include <sofa/simulation/task/ParallelForEach.h>
#include <sofa/simulation/task/TaskSchedulerUser.h>

namespace elasticity
{

template<class DataTypes, class ElementType>
class FEMForceField :
    public virtual sofa::core::behavior::ForceField<DataTypes>,
    public virtual sofa::core::behavior::TopologyAccessor,
    public virtual sofa::simulation::TaskSchedulerUser
{
public:
    SOFA_CLASS3(SOFA_TEMPLATE2(FEMForceField, DataTypes, ElementType),
        sofa::core::behavior::ForceField<DataTypes>,
        sofa::core::behavior::TopologyAccessor,
        sofa::simulation::TaskSchedulerUser);

private:
    using trait = elasticity::trait<DataTypes, ElementType>;
    using ElementForce = trait::ElementForce;

public:
    void init() override;

    void addForce(
        const sofa::core::MechanicalParams* mparams,
        sofa::DataVecDeriv_t<DataTypes>& f,
        const sofa::DataVecCoord_t<DataTypes>& x,
        const sofa::DataVecDeriv_t<DataTypes>& v) override;

    void addDForce(const sofa::core::MechanicalParams* mparams,
                   sofa::DataVecDeriv_t<DataTypes>& df,
                   const sofa::DataVecDeriv_t<DataTypes>& dx) override;

    void draw(const sofa::core::visual::VisualParams*) override;

    sofa::Data<ComputeStrategy> d_computeForceStrategy;
    sofa::Data<ComputeStrategy> d_computeForceDerivStrategy;

    sofa::Data<sofa::Real_t<DataTypes>> d_elementSpace;

protected:

    FEMForceField();

    void addElementForce(const sofa::core::MechanicalParams* mparams,
        sofa::type::vector<ElementForce>& f,
        const sofa::VecCoord_t<DataTypes>& x);

    virtual void beforeElementForce(const sofa::core::MechanicalParams* mparams,
        sofa::type::vector<ElementForce>& f,
        const sofa::VecCoord_t<DataTypes>& x) {}

    virtual void addElementForceRange(
        const sofa::simulation::Range<std::size_t>& range,
        const sofa::core::MechanicalParams* mparams,
        sofa::type::vector<ElementForce>& f,
        const sofa::VecCoord_t<DataTypes>& x) = 0;

    void addElementDForce(const sofa::core::MechanicalParams* mparams,
        sofa::type::vector<ElementForce>& df,
        const sofa::VecDeriv_t<DataTypes>& dx,
        sofa::Real_t<DataTypes> kFactor);

    virtual void addElementDForceRange(
        const sofa::simulation::Range<std::size_t>& range,
        const sofa::core::MechanicalParams* mparams,
        sofa::type::vector<ElementForce>& df,
        const sofa::VecDeriv_t<DataTypes>& dx,
        sofa::Real_t<DataTypes> kFactor) = 0;

    void dispatchElementForcesToNodes(
        const sofa::type::vector<typename trait::TopologyElement>& elements,
        sofa::VecDeriv_t<DataTypes>& nodeForces);

    sofa::type::vector<sofa::type::Vec<trait::NumberOfDofsInElement, sofa::Real_t<DataTypes>>> m_elementForce;
    sofa::type::vector<sofa::type::Vec<trait::NumberOfDofsInElement, sofa::Real_t<DataTypes>>> m_elementDForce;

    sofa::core::visual::DrawElementMesh<ElementType> m_drawMesh;
};


}  // namespace elasticity

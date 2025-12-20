#pragma once
#include <Elasticity/impl/ComputeStrategy.h>

namespace elasticity
{

template<class DataTypes, class ElementType>
class FEMForceField :
    public sofa::core::behavior::ForceField<DataTypes>,
    public virtual TopologyAccessor
{
public:
    SOFA_CLASS2(SOFA_TEMPLATE2(FEMForceField, DataTypes, ElementType),
        sofa::core::behavior::ForceField<DataTypes>,
        TopologyAccessor);

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

    sofa::Data<ComputeStrategy> d_computeForceStrategy;
    sofa::Data<ComputeStrategy> d_computeForceDerivStrategy;

protected:

    FEMForceField();

    void addElementForce(const sofa::core::MechanicalParams* mparams,
        sofa::type::vector<ElementForce>& f,
        const sofa::VecCoord_t<DataTypes>& x);

    void addElementDForce(const sofa::core::MechanicalParams* mparams,
        sofa::type::vector<ElementForce>& df,
        const sofa::VecDeriv_t<DataTypes>& dx,
        sofa::Real_t<DataTypes> kFactor);

    void dispatchElementForcesToNodes(
        const sofa::type::vector<typename trait::TopologyElement>& elements,
        sofa::VecDeriv_t<DataTypes>& nodeForces);

    sofa::type::vector<sofa::type::Vec<trait::NumberOfDofsInElement, sofa::Real_t<DataTypes>>> m_elementForce;
    sofa::type::vector<sofa::type::Vec<trait::NumberOfDofsInElement, sofa::Real_t<DataTypes>>> m_elementDForce;

    /**
     * The signature of the function that will be executed depending on the compute strategy of the
     * element forces
     */
    using ComputeElementForceFunction = std::function<void(sofa::type::vector<ElementForce>&, const sofa::VecCoord_t<DataTypes>&)>;

    /**
     * Stores the functions that will be executed depending on the compute strategy of the element
     * forces
     */
    std::unordered_map<std::string_view, ComputeElementForceFunction> m_computeElementForceMap;

    /**
     * The signature of the function that will be executed depending on the compute strategy of the
     * element forces derivatives
     */
    using ComputeElementForceDerivFunction = std::function<void(sofa::type::vector<ElementForce>&, const sofa::VecDeriv_t<DataTypes>&, sofa::Real_t<DataTypes>)>;

    /**
     * Stores the functions that will be executed depending on the compute strategy of the element
     * forces derivatives
     */
    std::unordered_map<std::string_view, ComputeElementForceDerivFunction> m_computeElementForceDerivMap;
};


}  // namespace elasticity

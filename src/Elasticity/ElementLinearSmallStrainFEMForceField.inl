#pragma once
#include <Elasticity/ElementLinearSmallStrainFEMForceField.h>
#include <Elasticity/MatrixTools.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/core/behavior/BaseLocalForceFieldMatrix.h>

namespace elasticity
{

template <class DataTypes, class ElementType>
ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::ElementLinearSmallStrainFEMForceField()
    : l_topology(initLink("topology", "Link to a topology containing elements"))
    , d_poissonRatio(
          initData(&d_poissonRatio, static_cast<Real>(0.45), "poissonRatio", "Poisson's ratio"))
    , d_youngModulus(
          initData(&d_youngModulus, static_cast<Real>(1e6), "youngModulus", "Young's modulus"))
{
    static std::string groupName = "Mechanical parameter";
    d_poissonRatio.setGroup(groupName);
    d_youngModulus.setGroup(groupName);
}

template <class DataTypes, class ElementType>
    void ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::init()
{
    Inherit1::init();

    if (!this->isComponentStateInvalid())
    {
        validateTopology();
    }

    if (!this->isComponentStateInvalid())
    {
        precomputeElementStiffness();
    }

    if (!this->isComponentStateInvalid())
    {
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Valid);
    }
}

template <class DataTypes, class ElementType>
void ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::validateTopology()
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

template <class DataTypes, class ElementType>
void ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::precomputeElementStiffness()
{
    const ElasticityTensor C = computeElasticityTensor();

    m_elementStiffness.clear();

    const auto& elements = FiniteElement::getElementSequence(*l_topology);
    m_elementStiffness.reserve(elements.size());

    auto restPositionAccessor = this->mstate->readRestPositions();
    for (const auto& element : elements)
    {
        //matrix where the i-th column is the i-th node coordinates in the element
        const sofa::type::Mat<spatial_dimensions, NumberOfNodesInElement, Real> X_element
            = nodesMatrix(element, restPositionAccessor.ref());

        ElementStiffness K;
        for (const auto& [quadraturePoint, weight] : FiniteElement::quadraturePoints())
        {
            //gradient of shape functions in the reference element evaluated at the quadrature point
            const sofa::type::Mat<NumberOfNodesInElement, ElementDimension, Real> dN_dq_ref = FiniteElement::gradientShapeFunctions(quadraturePoint);

            // jacobian of the mapping from the reference space to the physical space, evaluated at the quadrature point
            const sofa::type::Mat<spatial_dimensions, ElementDimension, Real> jacobian = X_element * dN_dq_ref;

            const auto detJ = elasticity::determinant(jacobian);
            const sofa::type::Mat<ElementDimension, spatial_dimensions, Real> J_inv = elasticity::inverse(jacobian);

            //gradient of the shape functions in the physical element evaluated at the quadrature point
            const sofa::type::Mat<NumberOfNodesInElement, spatial_dimensions, Real> dN_dq = dN_dq_ref * J_inv;

            const auto B = buildStrainDisplacement(dN_dq);

            K += (weight * detJ) * B.transposed() * C * B;
        }
        m_elementStiffness.push_back(K);
    }
}

template <class DataTypes, class ElementType>
void ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::addForce(
    const sofa::core::MechanicalParams* mparams,
    DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& v)
{
    SOFA_UNUSED(mparams);
    SOFA_UNUSED(v);

    auto forceAccessor = sofa::helper::getWriteOnlyAccessor(f);
    auto positionAccessor = sofa::helper::getReadAccessor(x);
    auto restPositionAccessor = this->mstate->readRestPositions();

    const auto& elements = FiniteElement::getElementSequence(*l_topology);

    Deriv nodeForce(sofa::type::NOINIT);

    auto elementStiffnessIt = m_elementStiffness.begin();
    for (const auto& element : elements)
    {
        std::array<Coord, NumberOfNodesInElement> elementNodesCoordinates, restElementNodesCoordinates;
        for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
        {
            elementNodesCoordinates[i] = positionAccessor[element[i]];
            restElementNodesCoordinates[i] = restPositionAccessor[element[i]];
        }

        const ElementDisplacement displacement = computeElementDisplacement(elementNodesCoordinates, restElementNodesCoordinates);

        const ElementStiffness& stiffnessMatrix = *elementStiffnessIt++;
        const sofa::type::Vec<NumberOfDofsInElement, Real> elementForce = stiffnessMatrix * displacement;

        for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
        {
            elementForce.getsub(i*spatial_dimensions, nodeForce);
            forceAccessor[element[i]] += -nodeForce;
        }
    }
}

template <class DataTypes, class ElementType>
void ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::addDForce(
    const sofa::core::MechanicalParams* mparams, DataVecDeriv& df, const DataVecDeriv& dx)
{
    auto dfAccessor = sofa::helper::getWriteAccessor(df);
    auto dxAccessor = sofa::helper::getReadAccessor(dx);
    dfAccessor.resize(dxAccessor.size());

    const Real kFactor = (Real)sofa::core::mechanicalparams::kFactorIncludingRayleighDamping(
        mparams, this->rayleighStiffness.getValue());

    const auto& elements = FiniteElement::getElementSequence(*l_topology);

    Deriv nodedForce(sofa::type::NOINIT);

    auto elementStiffnessIt = m_elementStiffness.begin();
    for (const auto& element : elements)
    {
        sofa::type::Vec<NumberOfDofsInElement, Real> element_dx;
        for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
        {
            for (sofa::Size j = 0; j < spatial_dimensions; ++j)
            {
                element_dx[i * spatial_dimensions + j] = dxAccessor[element[i]][j];
            }
        }

        const auto& stiffnessMatrix = *elementStiffnessIt++;
        const auto dForce = kFactor * stiffnessMatrix * element_dx;

        for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
        {
            dForce.getsub(i*spatial_dimensions, nodedForce);
            dfAccessor[element[i]] += -nodedForce;
        }
    }
}

template <class DataTypes, class ElementType>
void ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::buildStiffnessMatrix(
    sofa::core::behavior::StiffnessMatrix* matrix)
{
    sofa::type::Mat<spatial_dimensions, spatial_dimensions, Real> localMatrix(sofa::type::NOINIT);

    auto dfdx = matrix->getForceDerivativeIn(this->mstate)
                       .withRespectToPositionsIn(this->mstate);

    const auto& elements = FiniteElement::getElementSequence(*l_topology);
    auto elementStiffnessIt = m_elementStiffness.begin();
    for (const auto& element : elements)
    {
        const auto& stiffnessMatrix = *elementStiffnessIt++;

        for (sofa::Index n1 = 0; n1 < NumberOfNodesInElement; ++n1)
        {
            for (sofa::Index n2 = 0; n2 < NumberOfNodesInElement; ++n2)
            {
                stiffnessMatrix.getsub(spatial_dimensions * n1, spatial_dimensions * n2, localMatrix); //extract the submatrix corresponding to the coupling of nodes n1 and n2
                dfdx(element[n1] * spatial_dimensions, element[n2] * spatial_dimensions) += -localMatrix;
            }
        }
    }
}

template <class DataTypes, class ElementType>
SReal ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::getPotentialEnergy(
    const sofa::core::MechanicalParams*, const DataVecCoord& x) const
{
    return 0;
}

template <class DataTypes, class ElementType>
auto ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::computeElasticityTensor(
    Real youngModulus, Real poissonRatio)
-> ElasticityTensor
{
    static constexpr auto volumetricTensor = []()
    {
        sofa::type::Mat<NumberOfIndependentElements, NumberOfIndependentElements, Real> I_vol;
        for (sofa::Size i = 0; i < spatial_dimensions; ++i)
        {
            for (sofa::Size j = 0; j < spatial_dimensions; ++j)
            {
                I_vol[i][j] = 1;
            }
        }
        return I_vol;
    }();

    static constexpr auto deviatoricTensor = []()
    {
        sofa::type::Mat<NumberOfIndependentElements, NumberOfIndependentElements, Real> I_dev;
        for (sofa::Size i = 0; i < spatial_dimensions; ++i)
        {
            I_dev[i][i] = 1;
        }
        for (sofa::Size i = spatial_dimensions; i < NumberOfIndependentElements; ++i)
        {
            I_dev[i][i] = 0.5;
        }
        return I_dev;
    }();

    const auto mu = youngModulus / (2 * (1 + poissonRatio));
    const auto lambda = youngModulus * poissonRatio / ((1 + poissonRatio) * (1 - (spatial_dimensions - 1) * poissonRatio));

    return lambda * volumetricTensor + 2 * mu * deviatoricTensor;
}

template <class DataTypes, class ElementType>
auto ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::computeElasticityTensor() -> ElasticityTensor
{
    const auto E = d_youngModulus.getValue();
    const auto nu = d_poissonRatio.getValue();

    return computeElasticityTensor(E, nu);
}

template <class DataTypes, class ElementType>
typename ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::StrainDisplacement
ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::buildStrainDisplacement(
    const sofa::type::Mat<NumberOfNodesInElement, spatial_dimensions, Real> gradientShapeFunctions)
{
    StrainDisplacement B;
    for (sofa::Size ne = 0; ne < NumberOfNodesInElement; ++ne)
    {
        for (sofa::Size i = 0; i < spatial_dimensions; ++i)
        {
            B[i][ne * spatial_dimensions + i] = gradientShapeFunctions[ne][i];
        }

        auto row = spatial_dimensions;
        for (sofa::Size i = 0; i < spatial_dimensions; ++i)
        {
            for (sofa::Size j = i + 1; j < spatial_dimensions; ++j)
            {
                B[row][ne * spatial_dimensions + i] = gradientShapeFunctions[ne][j];
                B[row][ne * spatial_dimensions + j] = gradientShapeFunctions[ne][i];
                ++row;
            }
        }
    }

    return B;
}

template <class DataTypes, class ElementType>
auto ElementLinearSmallStrainFEMForceField<DataTypes, ElementType>::computeElementDisplacement(
    const std::array<Coord, NumberOfNodesInElement>& elementNodesCoordinates,
    const std::array<Coord, NumberOfNodesInElement>& restElementNodesCoordinates)
    -> ElementDisplacement
{
    ElementDisplacement displacement;
    for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
    {
        for (sofa::Size j = 0; j < spatial_dimensions; ++j)
        {
            displacement[i * spatial_dimensions + j] =
                elementNodesCoordinates[i][j] - restElementNodesCoordinates[i][j];
        }
    }
    return displacement;
}

}  // namespace elasticity

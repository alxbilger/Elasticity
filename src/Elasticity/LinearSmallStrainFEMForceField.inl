#pragma once
#include <Elasticity/LinearSmallStrainFEMForceField.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/core/behavior/BaseLocalForceFieldMatrix.h>
#include <Elasticity/FiniteElement.h>

namespace elasticity
{

template <class DataTypes, class ElementType>
LinearSmallStrainFEMForceField<DataTypes, ElementType>::LinearSmallStrainFEMForceField()
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
    void LinearSmallStrainFEMForceField<DataTypes, ElementType>::init()
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
void LinearSmallStrainFEMForceField<DataTypes, ElementType>::validateTopology()
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
void LinearSmallStrainFEMForceField<DataTypes, ElementType>::precomputeElementStiffness()
{
    const ElasticityTensor C = computeElasticityTensor();

    m_elementStiffness.clear();

    const auto& elements = FiniteElement::getElementSequence(*l_topology);
    m_elementStiffness.reserve(elements.size());

    auto restPositionAccessor = this->mstate->readRestPositions();
    for (const auto& element : elements)
    {
        std::array<Coord, NumberOfNodesInElement> restElementNodesCoordinates;
        for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
        {
            restElementNodesCoordinates[i] = restPositionAccessor[element[i]];
        }

        const Real volume = FiniteElement::volume(restElementNodesCoordinates);

        const auto B = computeStrainDisplacement(restElementNodesCoordinates);

        const auto K = volume * B.transposed() * C * B;
        m_elementStiffness.push_back(K);
    }
}

template <class DataTypes, class ElementType>
void LinearSmallStrainFEMForceField<DataTypes, ElementType>::addForce(
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
void LinearSmallStrainFEMForceField<DataTypes, ElementType>::addDForce(
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
void LinearSmallStrainFEMForceField<DataTypes, ElementType>::buildStiffnessMatrix(
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
SReal LinearSmallStrainFEMForceField<DataTypes, ElementType>::getPotentialEnergy(
    const sofa::core::MechanicalParams*, const DataVecCoord& x) const
{
    return 0;
}

template <class DataTypes, class ElementType>
auto LinearSmallStrainFEMForceField<DataTypes, ElementType>::computeElasticityTensor(
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
auto LinearSmallStrainFEMForceField<DataTypes, ElementType>::computeElasticityTensor() -> ElasticityTensor
{
    const auto E = d_youngModulus.getValue();
    const auto nu = d_poissonRatio.getValue();

    return computeElasticityTensor(E, nu);
}

template <class DataTypes, class ElementType>
auto LinearSmallStrainFEMForceField<DataTypes, ElementType>::computeShapeFunctions(
    const std::array<Coord, NumberOfNodesInElement>& elementNodesCoordinates)
-> std::array<ShapeFunction, NumberOfNodesInElement>
{
    sofa::type::Mat<NumberOfNodesInElement, NumberOfNodesInElement, Real> X;
    for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
    {
        X[i][0] = 1;
        for (sofa::Size j = 0; j < spatial_dimensions; ++j)
        {
            X[i][j+1] = elementNodesCoordinates[i][j];
        }
    }

    const auto invX = X.inverted();

    std::array<ShapeFunction, NumberOfNodesInElement> shapeFunctions;
    for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
    {
        shapeFunctions[i] = invX.col(i);
    }

    return shapeFunctions;
}

template <class DataTypes, class ElementType>
auto LinearSmallStrainFEMForceField<DataTypes, ElementType>::computeStrainDisplacement(
    const std::array<Coord, NumberOfNodesInElement>& elementNodesCoordinates) -> StrainDisplacement
{
    const std::array<ShapeFunction, NumberOfNodesInElement> shapeFunctions = computeShapeFunctions(elementNodesCoordinates);

    StrainDisplacement B;
    for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
    {
        B[0][i*spatial_dimensions] = shapeFunctions[i][1];
        B[1][i*spatial_dimensions + 1] = shapeFunctions[i][2];
        B[2][i*spatial_dimensions + 2] = shapeFunctions[i][3];

        B[3][i*spatial_dimensions + 1] = shapeFunctions[i][3];
        B[3][i*spatial_dimensions + 2] = shapeFunctions[i][2];

        B[4][i*spatial_dimensions + 0] = shapeFunctions[i][3];
        B[4][i*spatial_dimensions + 2] = shapeFunctions[i][1];

        B[5][i*spatial_dimensions + 0] = shapeFunctions[i][2];
        B[5][i*spatial_dimensions + 1] = shapeFunctions[i][1];
    }

    return B;
}

template <class DataTypes, class ElementType>
auto LinearSmallStrainFEMForceField<DataTypes, ElementType>::computeElementDisplacement(
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

#pragma once
#include <Elasticity/TET4LinearSmallStrainFEMForceField.h>
#include <sofa/core/ObjectFactory.h>

namespace elasticity
{

template <class DataTypes>
TET4LinearSmallStrainFEMForceField<DataTypes>::TET4LinearSmallStrainFEMForceField()
    : l_topology(initLink("topology", "Link to a topology containing tetrahedra"))
    , d_poissonRatio(
          initData(&d_poissonRatio, static_cast<Real>(0.45), "poissonRatio", "Poisson ratio"))
    , d_youngModulus(
          initData(&d_youngModulus, static_cast<Real>(1e6), "youngModulus", "Young's modulus"))
{
}

template <class DataTypes>
void TET4LinearSmallStrainFEMForceField<DataTypes>::init()
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

template <class DataTypes>
void TET4LinearSmallStrainFEMForceField<DataTypes>::validateTopology()
{
    if (l_topology.empty())
    {
        msg_info() << "link to Topology container should be set to ensure right behavior. First "
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

template <class DataTypes>
void TET4LinearSmallStrainFEMForceField<DataTypes>::precomputeElementStiffness()
{
    const auto C = computeElasticityTensor();

    m_elementStiffness.clear();

    const auto& tetrahedra = l_topology->getTetrahedra();
    m_elementStiffness.reserve(tetrahedra.size());

    auto restPositionAccessor = this->mstate->readRestPositions();
    for (const auto& tetrahedron : tetrahedra)
    {
        const auto [t0, t1, t2, t3] = tetrahedron.array();
        const auto volume =
            sofa::geometry::Tetrahedron::volume(restPositionAccessor[t0], restPositionAccessor[t1],
                                                restPositionAccessor[t2], restPositionAccessor[t3]);

        const auto B =
            computeStrainDisplacement({restPositionAccessor[t0], restPositionAccessor[t1],
                                       restPositionAccessor[t2], restPositionAccessor[t3]});

        const auto K = volume * B.transposed() * C * B;
        m_elementStiffness.push_back(K);
    }
}

template <class DataTypes>
void TET4LinearSmallStrainFEMForceField<DataTypes>::addForce(
    const sofa::core::MechanicalParams* mparams,
    DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& v)
{
    SOFA_UNUSED(mparams);
    SOFA_UNUSED(v);

    auto forceAccessor = sofa::helper::getWriteOnlyAccessor(f);
    auto positionAccessor = sofa::helper::getReadAccessor(x);
    auto restPositionAccessor = this->mstate->readRestPositions();

    const auto& tetrahedra = l_topology->getTetrahedra();

    auto elementStiffnessIt = m_elementStiffness.begin();
    for (const auto& tetrahedron : tetrahedra)
    {
        const auto [t0, t1, t2, t3] = tetrahedron.array();
        const std::array element_x { positionAccessor[t0], positionAccessor[t1], positionAccessor[t2], positionAccessor[t3]};
        const std::array element_X { restPositionAccessor[t0], restPositionAccessor[t1], restPositionAccessor[t2], restPositionAccessor[t3]};

        const auto displacement = computeElementDisplacement(element_x, element_X);

        const auto& K = *elementStiffnessIt++;
        const auto force = K * displacement;

        for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
        {
            for (sofa::Size j = 0; j < spatial_dimensions; ++j)
            {
                forceAccessor[tetrahedron[i]][j] += -force[i*3 + j];
            }
        }
    }
}

template <class DataTypes>
void TET4LinearSmallStrainFEMForceField<DataTypes>::addDForce(
    const sofa::core::MechanicalParams* mparams, DataVecDeriv& df, const DataVecDeriv& dx)
{
    auto dfAccessor = sofa::helper::getWriteAccessor(df);
    auto dxAccessor = sofa::helper::getReadAccessor(dx);
    dfAccessor.resize(dxAccessor.size());

    const Real kFactor = (Real)sofa::core::mechanicalparams::kFactorIncludingRayleighDamping(mparams, this->rayleighStiffness.getValue());

    const auto& tetrahedra = l_topology->getTetrahedra();

    auto elementStiffnessIt = m_elementStiffness.begin();
    for (const auto& tetrahedron : tetrahedra)
    {
        sofa::type::Vec<12, Real> element_dx;
        for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
        {
            for (sofa::Size j = 0; j < spatial_dimensions; ++j)
            {
                element_dx[i*3 + j] = dxAccessor[tetrahedron[i]][j];
            }
        }

        const auto& K = *elementStiffnessIt++;
        const auto dForce = kFactor * K * element_dx;

        for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
        {
            for (sofa::Size j = 0; j < spatial_dimensions; ++j)
            {
                dfAccessor[tetrahedron[i]][j] += -dForce[i*3 + j];
            }
        }
    }
}

template <class DataTypes>
SReal TET4LinearSmallStrainFEMForceField<DataTypes>::getPotentialEnergy(
    const sofa::core::MechanicalParams*, const DataVecCoord& x) const
{
    return 0;
}

template <class DataTypes>
auto TET4LinearSmallStrainFEMForceField<DataTypes>::computeElasticityTensor() -> ElasticityTensor
{
    const auto E = d_youngModulus.getValue();
    const auto nu = d_poissonRatio.getValue();

    ElasticityTensor C;

    const Real k = E / ((1 + nu) * (1 - 2 * nu));

    C(0, 0) = C(1, 1) = C(2, 2) = k * (1 - nu);
    C(3, 3) = C(4, 4) = C(5, 5) = k * (1 - 2 * nu) / 2;
    C(0, 1) = C(0, 2) = C(1, 0) = C(1, 2) = C(2, 0) = C(2, 1) = k * nu;

    return C;
}

template <class DataTypes>
auto TET4LinearSmallStrainFEMForceField<DataTypes>::computeStrainDisplacement(
    const std::array<Coord, NumberOfNodesInElement>& tetraNodesCoordinates) -> StrainDisplacement
{
    StrainDisplacement B;

    sofa::type::Mat<NumberOfNodesInElement, NumberOfNodesInElement, Real> X;
    for (int i = 0; i < NumberOfNodesInElement; ++i)
    {
        X[i][0] = tetraNodesCoordinates[i][0];
        X[i][1] = tetraNodesCoordinates[i][1];
        X[i][2] = tetraNodesCoordinates[i][2];
        X[i][3] = 1;
    }

    const auto invX = X.inverted();

    std::array<Coord, NumberOfNodesInElement> gradN;
    for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
    {
        for (sofa::Size j = 0; j < spatial_dimensions; ++j)
        {
            gradN[i][j] = invX[j+1][i];
        }
    }

    for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
    {
        B[0][i*spatial_dimensions] = gradN[i][0];
        B[1][i*spatial_dimensions + 1] = gradN[i][1];
        B[2][i*spatial_dimensions + 2] = gradN[i][2];
        B[3][i*spatial_dimensions] = gradN[i][1];
        B[3][i*spatial_dimensions + 1] = gradN[i][0];
        B[4][i*spatial_dimensions + 1] = gradN[i][2];
        B[4][i*spatial_dimensions + 2] = gradN[i][1];
        B[5][i*spatial_dimensions] = gradN[i][2];
        B[5][i*spatial_dimensions + 2] = gradN[i][0];
    }

    return B;
}

template <class DataTypes>
auto TET4LinearSmallStrainFEMForceField<DataTypes>::computeElementDisplacement(
    const std::array<Coord, NumberOfNodesInElement>& tetraNodesCoordinates,
    const std::array<Coord, NumberOfNodesInElement>& restTetraNodesCoordinates) -> ElementDisplacement
{
    ElementDisplacement displacement;
    for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
    {
        for (sofa::Size j = 0; j < spatial_dimensions; ++j)
        {
            displacement[i*spatial_dimensions + j] = tetraNodesCoordinates[i][j] - restTetraNodesCoordinates[i][j];
        }
    }
    return displacement;
}

}  // namespace elasticity

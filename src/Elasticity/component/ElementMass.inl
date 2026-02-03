#pragma once
#include <Elasticity/component/ElementMass.h>
#include <Elasticity/impl/VectorTools.h>
#include <Elasticity/impl/MatrixTools.h>

namespace elasticity
{

template <class DataTypes, class ElementType>
ElementMass<DataTypes, ElementType>::ElementMass()
    : d_nodalDensity(initData(&d_nodalDensity, {defaultNodalDensity}, "nodalDensity", "Scalar density assigned to each node. The density is interpolated between nodes using the element shape functions."))
{
}

template <class DataTypes, class ElementType>
void ElementMass<DataTypes, ElementType>::init()
{
    sofa::core::behavior::Mass<DataTypes>::init();

    if (!this->isComponentStateInvalid())
    {
        sofa::core::behavior::TopologyAccessor::init();
    }

    if (!this->isComponentStateInvalid() && this->mstate)
    {
        this->resizeNodalDensity(this->mstate->getSize());
    }

    // (Re)build cache if possible
    if (!this->isComponentStateInvalid())
    {
        updateMassCacheIfNeeded();
    }

    if (!this->isComponentStateInvalid())
    {
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Valid);
    }
}

template <class DataTypes, class ElementType>
void ElementMass<DataTypes, ElementType>::resizeNodalDensity(const std::size_t size)
{
    sofa::helper::WriteAccessor nodalDensity = sofa::helper::getWriteAccessor(d_nodalDensity);

    if (nodalDensity.size() < size)
    {
        if (nodalDensity.empty())
        {
            nodalDensity.resize(size, defaultNodalDensity);
        }
        else
        {
            nodalDensity.resize(size, nodalDensity.back());
        }
    }
}

template <class DataTypes, class ElementType>
sofa::Deriv_t<DataTypes> ElementMass<DataTypes, ElementType>::getGravity() const
{
    const auto& nodeGravity = this->getContext()->getGravity();

    sofa::Deriv_t<DataTypes> gravity;
    DataTypes::set(gravity, nodeGravity[0], nodeGravity[1], nodeGravity[2]);

    return gravity;
}

template <class DataTypes, class ElementType>
bool ElementMass<DataTypes, ElementType>::isDiagonal() const
{
    return false;
}

template <class DataTypes, class ElementType>
void ElementMass<DataTypes, ElementType>::updateMassCacheIfNeeded()
{
    if (!this->l_topology || !this->mstate)
        return;

    const auto& elements = FiniteElement::getElementSequence(*this->l_topology);
    if (m_elementQuadratureMass.size() == elements.size())
        return;

    m_elementQuadratureMass.clear();
    m_elementQuadratureMass.resize(elements.size());

    const auto restPositionsAccessor = this->mstate->readRestPositions();

    std::size_t elementId = 0;
    for (const auto& element : elements)
    {
        const std::array<Coord, NumberOfNodesInElement> elementNodesRestCoordinates =
            extractNodesVectorFromGlobalVector(element, restPositionsAccessor.ref());

        auto& elementData = m_elementQuadratureMass[elementId];

        std::size_t qIndex = 0;
        for (const auto& [quadraturePoint, weight] : FiniteElement::quadraturePoints())
        {
            const auto N = FiniteElement::shapeFunctions(quadraturePoint);
            const auto dN_dq_ref = FiniteElement::gradientShapeFunctions(quadraturePoint);

            sofa::type::Mat<spatial_dimensions, FiniteElement::TopologicalDimension, Real> jacobian;
            for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
            {
                jacobian += sofa::type::dyad(elementNodesRestCoordinates[i], dN_dq_ref[i]);
            }

            const auto detJ = elasticity::absGeneralizedDeterminant(jacobian);
            const Real coeff = static_cast<Real>(weight) * static_cast<Real>(detJ);

            ScalarMassMatrix& Bq = elementData[qIndex];
            Bq.clear();

            for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
            {
                for (sofa::Size j = 0; j < NumberOfNodesInElement; ++j)
                {
                    Bq(i, j) = coeff * N[i] * N[j];
                }
            }

            ++qIndex;
        }

        ++elementId;
    }
}

template <class DataTypes, class ElementType>
void ElementMass<DataTypes, ElementType>::addMDx(const sofa::core::MechanicalParams* mparams,
                                                sofa::DataVecDeriv_t<DataTypes>& vres,
                                                const sofa::DataVecDeriv_t<DataTypes>& vdx,
                                                SReal factor)
{
    SOFA_UNUSED(mparams);

    if (!this->l_topology)
        return;

    updateMassCacheIfNeeded();

    const sofa::helper::ReadAccessor nodalDensity = sofa::helper::getReadAccessor(d_nodalDensity);
    auto res = sofa::helper::getWriteAccessor(vres);
    auto dx = sofa::helper::getReadAccessor(vdx);
    resizeNodalDensity(dx.size());

    if (res.size() < dx.size())
        res.resize(dx.size());

    const auto& elements = FiniteElement::getElementSequence(*this->l_topology);
    if (m_elementQuadratureMass.size() != elements.size())
        return;

    const Real fact = static_cast<Real>(factor);

    for (std::size_t elementId = 0; elementId < elements.size(); ++elementId)
    {
        const auto& element = elements[elementId];

        const std::array<Real, NumberOfNodesInElement> elementNodesDensity =
            extractNodesVectorFromGlobalVector(element, nodalDensity.ref());

        const auto& elementData = m_elementQuadratureMass[elementId];

        std::size_t qIndex = 0;
        for (const auto& [quadraturePoint, weight] : FiniteElement::quadraturePoints())
        {
            SOFA_UNUSED(weight);

            const auto N = FiniteElement::shapeFunctions(quadraturePoint);

            Real rho_q = static_cast<Real>(0);
            for (sofa::Size k = 0; k < NumberOfNodesInElement; ++k)
            {
                rho_q += N[k] * elementNodesDensity[k];
            }

            const ScalarMassMatrix& Bq = elementData[qIndex];

            for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
            {
                Deriv acc{};
                for (sofa::Size j = 0; j < NumberOfNodesInElement; ++j)
                {
                    acc += dx[element[j]] * (rho_q * Bq(i, j));
                }
                res[element[i]] += acc * fact;
            }

            ++qIndex;
        }
    }
}

template <class DataTypes, class ElementType>
void ElementMass<DataTypes, ElementType>::addMToMatrix(sofa::linearalgebra::BaseMatrix * matrix, SReal mFact, unsigned int &offset)
{
    if (!matrix || !this->l_topology || !this->mstate)
        return;

    updateMassCacheIfNeeded();

    const sofa::helper::ReadAccessor nodalDensity = sofa::helper::getReadAccessor(d_nodalDensity);
    resizeNodalDensity(this->mstate->getSize());

    const auto& elements = FiniteElement::getElementSequence(*this->l_topology);
    if (m_elementQuadratureMass.size() != elements.size())
        return;

    // Assemble: for each element, for each quadrature point:
    // add (mFact * rho_q * Bq(i,j)) on each spatial component (block diagonal)
    for (std::size_t elementId = 0; elementId < elements.size(); ++elementId)
    {
        const auto& element = elements[elementId];

        const std::array<Real, NumberOfNodesInElement> elementNodesDensity =
            extractNodesVectorFromGlobalVector(element, nodalDensity.ref());

        const auto& elementData = m_elementQuadratureMass[elementId];

        std::size_t qIndex = 0;
        for (const auto& [quadraturePoint, weight] : FiniteElement::quadraturePoints())
        {
            SOFA_UNUSED(weight);

            const auto N = FiniteElement::shapeFunctions(quadraturePoint);

            Real rho_q = static_cast<Real>(0);
            for (sofa::Size k = 0; k < NumberOfNodesInElement; ++k)
            {
                rho_q += N[k] * elementNodesDensity[k];
            }

            const ScalarMassMatrix& Bq = elementData[qIndex];

            for (sofa::Size i = 0; i < NumberOfNodesInElement; ++i)
            {
                const sofa::Index gi = element[i];
                for (sofa::Size j = 0; j < NumberOfNodesInElement; ++j)
                {
                    const sofa::Index gj = element[j];

                    const Real mij = mFact * rho_q * Bq(i, j);
                    if (mij == static_cast<Real>(0))
                        continue;

                    for (sofa::Size d = 0; d < spatial_dimensions; ++d)
                    {
                        const auto row = static_cast<sofa::linearalgebra::BaseMatrix::Index>(offset + gi * spatial_dimensions + d);
                        const auto col = static_cast<sofa::linearalgebra::BaseMatrix::Index>(offset + gj * spatial_dimensions + d);
                        matrix->add(row, col, static_cast<SReal>(mij));
                    }
                }
            }

            ++qIndex;
        }
    }
}

template <class DataTypes, class ElementType>
void ElementMass<DataTypes, ElementType>::addForce(const sofa::core::MechanicalParams* mparams,
                                                   sofa::DataVecDeriv_t<DataTypes>& f,
                                                   const sofa::DataVecCoord_t<DataTypes>& x,
                                                   const sofa::DataVecDeriv_t<DataTypes>& v)
{
    SOFA_UNUSED(mparams);
    SOFA_UNUSED(v);

    if (!l_topology)
        return;

    const sofa::helper::ReadAccessor nodalDensity = sofa::helper::getReadAccessor(d_nodalDensity);

    auto positionAccessor = sofa::helper::getReadAccessor(x);
    auto forceAccessor = sofa::helper::getWriteAccessor(f);

    //make sure there is a density associated to each node
    resizeNodalDensity(positionAccessor.size());

    const auto gravity = this->getGravity();

    const auto& elements = FiniteElement::getElementSequence(*this->l_topology);

    for (std::size_t elementId = 0; elementId < elements.size(); ++elementId)
    {
        const auto& element = elements[elementId];

        const std::array<sofa::Coord_t<DataTypes>, trait::NumberOfNodesInElement> elementNodesCoordinates =
            extractNodesVectorFromGlobalVector(element, positionAccessor.ref());

        const std::array<sofa::Real_t<DataTypes>, trait::NumberOfNodesInElement> elementNodesDensity =
            extractNodesVectorFromGlobalVector(element, nodalDensity.ref());

        // Integrate: f_i += ∫_Ω rho(x) * N_i(x) * g dΩ
        // with rho(x) interpolated from nodal densities: rho(x) = Σ_j N_j(x) * rho_j
        for (const auto& [quadraturePoint, weight] : FiniteElement::quadraturePoints())
        {
            // shape functions evaluated at the quadrature point
            const auto N = FiniteElement::shapeFunctions(quadraturePoint);

            // gradient of shape functions in the reference element evaluated at the quadrature point
            const auto dN_dq_ref = FiniteElement::gradientShapeFunctions(quadraturePoint);

            // jacobian of the mapping from the reference space to the physical space
            sofa::type::Mat<DataTypes::spatial_dimensions, FiniteElement::TopologicalDimension, sofa::Real_t<DataTypes>> jacobian;
            for (sofa::Size i = 0; i < trait::NumberOfNodesInElement; ++i)
            {
                jacobian += sofa::type::dyad(elementNodesCoordinates[i], dN_dq_ref[i]);
            }

            const auto detJ = elasticity::absGeneralizedDeterminant(jacobian);
            const auto dV = weight * detJ;

            // density at quadrature point (interpolated)
            sofa::Real_t<DataTypes> rho_q = static_cast<sofa::Real_t<DataTypes>>(0);
            for (sofa::Size j = 0; j < trait::NumberOfNodesInElement; ++j)
            {
                rho_q += N[j] * elementNodesDensity[j];
            }

            // distribute body force to element nodes
            for (sofa::Size i = 0; i < trait::NumberOfNodesInElement; ++i)
            {
                forceAccessor[element[i]] += (dV * (rho_q * N[i])) * gravity;
            }
        }
    }
}

}  // namespace elasticity

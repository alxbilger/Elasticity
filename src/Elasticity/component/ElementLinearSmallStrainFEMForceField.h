#pragma once

#include <Elasticity/component/LinearMechanicalParametersComponent.h>
#include <Elasticity/component/TopologyAccessor.h>
#include <Elasticity/config.h>
#include <Elasticity/impl/ElementStiffnessMatrix.h>
#include <Elasticity/impl/trait.h>
#include <sofa/core/behavior/ForceField.h>
#include <sofa/helper/OptionsGroup.h>

#if !defined(ELASTICITY_COMPONENT_ELEMENT_LINEAR_SMALL_STRAIN_FEM_FORCE_FIELD_CPP)
#include <Elasticity/finiteelement/FiniteElement[all].h>
#endif

namespace elasticity
{

template <class DataTypes, class ElementType>
struct ComputeElementForceStrategy;

template <class DataTypes, class ElementType>
class ElementLinearSmallStrainFEMForceField :
    public TopologyAccessor,
    public LinearMechanicalParametersComponent<DataTypes>,
    public sofa::core::behavior::ForceField<DataTypes>
{
public:
    SOFA_CLASS3(
        SOFA_TEMPLATE2(ElementLinearSmallStrainFEMForceField, DataTypes, ElementType),
            TopologyAccessor,
            LinearMechanicalParametersComponent<DataTypes>,
            sofa::core::behavior::ForceField<DataTypes>);

    /**
     * The purpose of this function is to register the name of this class according to the provided
     * pattern.
     *
     * Example: ElementLinearSmallStrainFEMForceField<Vec3Types, sofa::geometry::Edge> will produce
     * the class name "EdgeLinearSmallStrainFEMForceField".
     */
    static const std::string GetCustomClassName()
    {
        return std::string(sofa::geometry::elementTypeToString(ElementType::Element_type)) +
               "LinearSmallStrainFEMForceField";
    }

    static const std::string GetCustomTemplateName() { return DataTypes::Name(); }

    sofa::Data<sofa::helper::OptionsGroup> d_computeStrategy;

private:
    using trait = elasticity::trait<DataTypes, ElementType>;
    using ElementStiffness = typename trait::ElementStiffness;
    using ElasticityTensor = typename trait::ElasticityTensor;
    using ElementDisplacement = typename trait::ElementDisplacement;
    using StrainDisplacement = typename trait::StrainDisplacement;

public:
    void init() override;

    void addForce(const sofa::core::MechanicalParams* mparams, sofa::DataVecDeriv_t<DataTypes>& f,
              const sofa::DataVecCoord_t<DataTypes>& x, const sofa::DataVecDeriv_t<DataTypes>& v) override;

    void addDForce(const sofa::core::MechanicalParams* mparams, sofa::DataVecDeriv_t<DataTypes>& df,
                   const sofa::DataVecDeriv_t<DataTypes>& dx) override;

    void buildStiffnessMatrix(sofa::core::behavior::StiffnessMatrix* matrix) override;

    SReal getPotentialEnergy(const sofa::core::MechanicalParams*, const sofa::DataVecCoord_t<DataTypes>& x) const override;

    using sofa::core::behavior::ForceField<DataTypes>::addKToMatrix;
    // almost deprecated, but here for compatibility with unit tests
    void addKToMatrix(sofa::linearalgebra::BaseMatrix* matrix, SReal kFact, unsigned& offset) override;

protected:

    ElementLinearSmallStrainFEMForceField();

    /**
     * With linear small strain, the element stiffness matrix is constant, so it can be precomputed.
     */
    void precomputeElementStiffness();

    /**
     * List of precomputed element stiffness matrices
     */
    sofa::type::vector<ElementStiffness> m_elementStiffness;

    ElasticityTensor m_elasticityTensor;

    sofa::type::vector<std::array<StrainDisplacement, trait::NumberOfNodesInElement>> m_strainDisplacement;

    template <typename TCoord>
    static ElementDisplacement computeElementDisplacement(
        const std::array<TCoord, trait::NumberOfNodesInElement>& elementNodesCoordinates,
        const std::array<TCoord, trait::NumberOfNodesInElement>& restElementNodesCoordinates)
    {
        ElementDisplacement displacement(sofa::type::NOINIT);
        for (sofa::Size i = 0; i < trait::NumberOfNodesInElement; ++i)
        {
            for (sofa::Size j = 0; j < trait::spatial_dimensions; ++j)
            {
                displacement[i * trait::spatial_dimensions + j] = elementNodesCoordinates[i][j] - restElementNodesCoordinates[i][j];
            }
        }
        return displacement;
    }

    void selectStrategy();

    std::unique_ptr<ComputeElementForceStrategy<DataTypes, ElementType>> m_computeElementForceStrategy;

    sofa::type::vector<sofa::type::Vec<trait::NumberOfDofsInElement, sofa::Real_t<DataTypes>>> m_elementForce;
};

template <class DataTypes, class ElementType>
struct ComputeElementForceStrategy
{
    using trait = elasticity::trait<DataTypes, ElementType>;
    using TopologyElement = typename trait::TopologyElement;
    using ElementStiffness = typename trait::ElementStiffness;
    using Real = sofa::Real_t<DataTypes>;

    virtual ~ComputeElementForceStrategy() = default;
    virtual void compute(
        const sofa::type::vector<TopologyElement>& elements,
        const sofa::VecCoord_t<DataTypes>& position,
        const sofa::VecCoord_t<DataTypes>& restPosition,
        sofa::type::vector<sofa::type::Vec<trait::NumberOfDofsInElement, Real>>& elementForces) = 0;

    virtual void setElementStiffnessMatrices(const sofa::type::vector<ElementStiffness>& m_elementStiffness) = 0;
};

#if !defined(ELASTICITY_COMPONENT_ELEMENT_LINEAR_SMALL_STRAIN_FEM_FORCE_FIELD_CPP)
extern template class ELASTICITY_API ElementLinearSmallStrainFEMForceField<sofa::defaulttype::Vec1Types, sofa::geometry::Edge>;
extern template class ELASTICITY_API ElementLinearSmallStrainFEMForceField<sofa::defaulttype::Vec2Types, sofa::geometry::Edge>;
extern template class ELASTICITY_API ElementLinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Edge>;
extern template class ELASTICITY_API ElementLinearSmallStrainFEMForceField<sofa::defaulttype::Vec2Types, sofa::geometry::Triangle>;
extern template class ELASTICITY_API ElementLinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Triangle>;
extern template class ELASTICITY_API ElementLinearSmallStrainFEMForceField<sofa::defaulttype::Vec2Types, sofa::geometry::Quad>;
extern template class ELASTICITY_API ElementLinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Quad>;
extern template class ELASTICITY_API ElementLinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Tetrahedron>;
extern template class ELASTICITY_API ElementLinearSmallStrainFEMForceField<sofa::defaulttype::Vec3Types, sofa::geometry::Hexahedron>;
#endif

}  // namespace elasticity

#pragma once

#include <sofa/core/behavior/BaseLocalForceFieldMatrix.h>

namespace elasticity
{

/**
 * Base class for FEM. The class is designed to be called in a ForceField.
 */
template <class DataTypes>
class BaseFEM
{
    using VecCoord = sofa::VecCoord_t<DataTypes>;
    using VecDeriv = sofa::VecDeriv_t<DataTypes>;
    using Real = sofa::Real_t<DataTypes>;

public:
    virtual ~BaseFEM() = default;

    virtual void addForce(VecDeriv& force, const VecCoord& position, const VecCoord& restPosition) = 0;
    virtual void addDForce(VecDeriv& df, const VecDeriv& dx, Real kFactor) const = 0;
    virtual void buildStiffnessMatrix(sofa::core::behavior::StiffnessMatrix::Derivative& dfdx) const = 0;
};

}

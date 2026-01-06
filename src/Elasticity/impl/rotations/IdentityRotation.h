#pragma once

namespace elasticity
{

struct IdentityRotation
{
    template<class RotationMatrix, class NotUsed1, class NotUsed2>
    void computeRotation(RotationMatrix& rotationMatrix, const NotUsed1&, const NotUsed2&)
    {
        rotationMatrix.identity();
    }

    static constexpr sofa::helper::Item getItem()
    {
        return {"identity", "Identity rotation. Equivalent to the linear small strain FEM."};
    }
};

}

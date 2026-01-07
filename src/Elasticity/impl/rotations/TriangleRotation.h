#pragma once

#include <sofa/helper/SelectableItem.h>
#include <sofa/type/Mat.h>

namespace elasticity
{

template <class DataTypes>
struct TriangleRotation
{
    using RotationMatrix = sofa::type::Mat<DataTypes::spatial_dimensions, DataTypes::spatial_dimensions, sofa::Real_t<DataTypes>>;

    template<sofa::Size NumberOfNodesInElement>
    void computeRotation(RotationMatrix& rotationMatrix, const RotationMatrix& initialRotationMatrix,
        const std::array<sofa::Coord_t<DataTypes>, NumberOfNodesInElement>& nodesPosition,
        const std::array<sofa::Coord_t<DataTypes>, NumberOfNodesInElement>& nodesRestPosition)
    {
        SOFA_UNUSED(nodesRestPosition);

        RotationMatrix currentRotation(sofa::type::NOINIT);
        computeRotationFrom3Points(currentRotation, {nodesPosition[0], nodesPosition[1], nodesPosition[2]});

        rotationMatrix = currentRotation.multTranspose(initialRotationMatrix);
    }

    static constexpr sofa::helper::Item getItem()
    {
        return {"triangle", "Compute the rotation based on the Gram-Schmidt frame alignment"};
    }

private:

    void computeRotationFrom3Points(RotationMatrix& rotationMatrix,
        const std::array<sofa::Coord_t<DataTypes>, 3>& nodesPosition)
    {
        using Coord = sofa::Coord_t<DataTypes>;

        const Coord edgex = (nodesPosition[1] - nodesPosition[0]).normalized();
        Coord edgey = nodesPosition[2] - nodesPosition[0];
        const Coord edgez = cross( edgex, edgey ).normalized();
        edgey = cross( edgez, edgex ); //edgey is unit vector because edgez and edgex are orthogonal unit vectors

        rotationMatrix(0,0) = edgex[0];
        rotationMatrix(0,1) = edgex[1];
        rotationMatrix(0,2) = edgex[2];

        rotationMatrix(1,0) = edgey[0];
        rotationMatrix(1,1) = edgey[1];
        rotationMatrix(1,2) = edgey[2];

        rotationMatrix(2,0) = edgez[0];
        rotationMatrix(2,1) = edgez[1];
        rotationMatrix(2,2) = edgez[2];
    }
};

}

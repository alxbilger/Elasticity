#pragma once

#include <sofa/core/trait/DataTypes.h>
#include <sofa/helper/SelectableItem.h>
#include <sofa/type/Mat.h>

namespace elasticity
{

template <class DataTypes>
struct HexaRotation
{
    using RotationMatrix = sofa::type::Mat<DataTypes::spatial_dimensions, DataTypes::spatial_dimensions, sofa::Real_t<DataTypes>>;

    template<sofa::Size NumberOfNodesInElement>
    void computeRotation(RotationMatrix& rotationMatrix, const RotationMatrix& initialRotationMatrix,
        const std::array<sofa::Coord_t<DataTypes>, NumberOfNodesInElement>& nodesPosition,
        const std::array<sofa::Coord_t<DataTypes>, NumberOfNodesInElement>& nodesRestPosition)
    {
        SOFA_UNUSED(nodesRestPosition);

        RotationMatrix currentRotation(sofa::type::NOINIT);
        computeRotationFromHexa(currentRotation, nodesPosition);

        rotationMatrix = currentRotation.multTranspose(initialRotationMatrix);
    }

    static constexpr sofa::helper::Item getItem()
    {
        return {"hexa", "Compute the rotation based on two average edges in the hexahedron"};
    }

private:

    template<sofa::Size NumberOfNodesInElement>
    void computeRotationFromHexa(RotationMatrix& rotationMatrix,
        const std::array<sofa::Coord_t<DataTypes>, NumberOfNodesInElement>& nodes)
    {
        using Coord = sofa::Coord_t<DataTypes>;

        Coord edgex = (nodes[1]-nodes[0] + nodes[2]-nodes[3] + nodes[5]-nodes[4] + nodes[6]-nodes[7])*.25;
        Coord edgey = (nodes[3]-nodes[0] + nodes[2]-nodes[1] + nodes[7]-nodes[4] + nodes[6]-nodes[5])*.25;

        edgex.normalize();

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

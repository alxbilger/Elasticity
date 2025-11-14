#include <Elasticity/CUDA/config.h>

#include <Elasticity/component/ElementLinearSmallStrainFEMForceField.inl>
#include <sofa/core/ObjectFactory.h>
#include <sofa/gpu/cuda/CudaTypes.h>

namespace elasticity
{

template class ELASTICITY_CUDA_API ElementLinearSmallStrainFEMForceField<sofa::gpu::cuda::CudaVec1fTypes, sofa::geometry::Edge>;
template class ELASTICITY_CUDA_API ElementLinearSmallStrainFEMForceField<sofa::gpu::cuda::CudaVec2fTypes, sofa::geometry::Edge>;
template class ELASTICITY_CUDA_API ElementLinearSmallStrainFEMForceField<sofa::gpu::cuda::CudaVec3fTypes, sofa::geometry::Edge>;
template class ELASTICITY_CUDA_API ElementLinearSmallStrainFEMForceField<sofa::gpu::cuda::CudaVec2fTypes, sofa::geometry::Triangle>;
template class ELASTICITY_CUDA_API ElementLinearSmallStrainFEMForceField<sofa::gpu::cuda::CudaVec3fTypes, sofa::geometry::Triangle>;
template class ELASTICITY_CUDA_API ElementLinearSmallStrainFEMForceField<sofa::gpu::cuda::CudaVec2fTypes, sofa::geometry::Quad>;
template class ELASTICITY_CUDA_API ElementLinearSmallStrainFEMForceField<sofa::gpu::cuda::CudaVec3fTypes, sofa::geometry::Quad>;
template class ELASTICITY_CUDA_API ElementLinearSmallStrainFEMForceField<sofa::gpu::cuda::CudaVec3fTypes, sofa::geometry::Tetrahedron>;
template class ELASTICITY_CUDA_API ElementLinearSmallStrainFEMForceField<sofa::gpu::cuda::CudaVec3fTypes, sofa::geometry::Hexahedron>;

namespace cuda
{

void registerCudaElementLinearSmallStrainFEMForceField(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData("Hooke's law on linear beams assuming small strain")
            .add< ElementLinearSmallStrainFEMForceField<sofa::gpu::cuda::CudaVec1fTypes, sofa::geometry::Edge> >()
            .add< ElementLinearSmallStrainFEMForceField<sofa::gpu::cuda::CudaVec2fTypes, sofa::geometry::Edge> >()
            .add< ElementLinearSmallStrainFEMForceField<sofa::gpu::cuda::CudaVec3fTypes, sofa::geometry::Edge> >());

    factory->registerObjects(sofa::core::ObjectRegistrationData("Hooke's law on linear triangles assuming small strain")
        .add< ElementLinearSmallStrainFEMForceField<sofa::gpu::cuda::CudaVec2fTypes, sofa::geometry::Triangle> >()
        .add< ElementLinearSmallStrainFEMForceField<sofa::gpu::cuda::CudaVec3fTypes, sofa::geometry::Triangle> >());

    factory->registerObjects(sofa::core::ObjectRegistrationData("Hooke's law on linear quads assuming small strain")
        .add< ElementLinearSmallStrainFEMForceField<sofa::gpu::cuda::CudaVec2fTypes, sofa::geometry::Quad> >()
        .add< ElementLinearSmallStrainFEMForceField<sofa::gpu::cuda::CudaVec3fTypes, sofa::geometry::Quad> >());

    factory->registerObjects(sofa::core::ObjectRegistrationData("Hooke's law on linear tetrahedra assuming small strain")
        .add< ElementLinearSmallStrainFEMForceField<sofa::gpu::cuda::CudaVec3fTypes, sofa::geometry::Tetrahedron> >());

    factory->registerObjects(sofa::core::ObjectRegistrationData("Hooke's law on linear hexahedra assuming small strain")
        .add< ElementLinearSmallStrainFEMForceField<sofa::gpu::cuda::CudaVec3fTypes, sofa::geometry::Hexahedron> >());
}

}


}

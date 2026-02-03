#include <Elasticity/CUDA/config.h>

#include <Elasticity/component/ElementCorotationalFEMForceField.inl>
#include <sofa/core/ObjectFactory.h>
#include <sofa/gpu/cuda/CudaTypes.h>

namespace elasticity
{

// template class ELASTICITY_CUDA_API ElementCorotationalFEMForceField<sofa::gpu::cuda::CudaVec1fTypes, sofa::geometry::Edge>;
template class ELASTICITY_CUDA_API ElementCorotationalFEMForceField<sofa::gpu::cuda::CudaVec2fTypes, sofa::geometry::Edge>;
template class ELASTICITY_CUDA_API ElementCorotationalFEMForceField<sofa::gpu::cuda::CudaVec3fTypes, sofa::geometry::Edge>;
template class ELASTICITY_CUDA_API ElementCorotationalFEMForceField<sofa::gpu::cuda::CudaVec2fTypes, sofa::geometry::Triangle>;
template class ELASTICITY_CUDA_API ElementCorotationalFEMForceField<sofa::gpu::cuda::CudaVec3fTypes, sofa::geometry::Triangle>;
template class ELASTICITY_CUDA_API ElementCorotationalFEMForceField<sofa::gpu::cuda::CudaVec2fTypes, sofa::geometry::Quad>;
template class ELASTICITY_CUDA_API ElementCorotationalFEMForceField<sofa::gpu::cuda::CudaVec3fTypes, sofa::geometry::Quad>;
template class ELASTICITY_CUDA_API ElementCorotationalFEMForceField<sofa::gpu::cuda::CudaVec3fTypes, sofa::geometry::Tetrahedron>;
template class ELASTICITY_CUDA_API ElementCorotationalFEMForceField<sofa::gpu::cuda::CudaVec3fTypes, sofa::geometry::Hexahedron>;


namespace cuda
{

void registerCudaElementCorotationalFEMForceField(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData("Hooke's law on linear beams using the corotational approach")
//     .add< ElementCorotationalFEMForceField<sofa::gpu::cuda::CudaVec1fTypes, sofa::geometry::Edge> >()
        .add< ElementCorotationalFEMForceField<sofa::gpu::cuda::CudaVec2fTypes, sofa::geometry::Edge> >()
        .add< ElementCorotationalFEMForceField<sofa::gpu::cuda::CudaVec3fTypes, sofa::geometry::Edge> >());

    factory->registerObjects(sofa::core::ObjectRegistrationData("Hooke's law on linear triangles using the corotational approach")
        .add< ElementCorotationalFEMForceField<sofa::gpu::cuda::CudaVec2fTypes, sofa::geometry::Triangle> >()
        .add< ElementCorotationalFEMForceField<sofa::gpu::cuda::CudaVec3fTypes, sofa::geometry::Triangle> >());

    factory->registerObjects(sofa::core::ObjectRegistrationData("Hooke's law on linear quads using the corotational approach")
        .add< ElementCorotationalFEMForceField<sofa::gpu::cuda::CudaVec2fTypes, sofa::geometry::Quad> >()
        .add< ElementCorotationalFEMForceField<sofa::gpu::cuda::CudaVec3fTypes, sofa::geometry::Quad> >());

    factory->registerObjects(sofa::core::ObjectRegistrationData("Hooke's law on linear tetrahedra using the corotational approach")
        .add< ElementCorotationalFEMForceField<sofa::gpu::cuda::CudaVec3fTypes, sofa::geometry::Tetrahedron> >());

    factory->registerObjects(sofa::core::ObjectRegistrationData("Hooke's law on linear hexahedra using the corotational approach")
        .add< ElementCorotationalFEMForceField<sofa::gpu::cuda::CudaVec3fTypes, sofa::geometry::Hexahedron> >());
}

}


}

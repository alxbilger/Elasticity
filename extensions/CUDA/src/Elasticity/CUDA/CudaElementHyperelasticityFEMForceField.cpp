#include <Elasticity/CUDA/config.h>

#include <Elasticity/component/ElementHyperelasticityFEMForceField.inl>
#include <sofa/core/ObjectFactory.h>
#include <sofa/gpu/cuda/CudaTypes.h>

namespace elasticity
{

template class ELASTICITY_CUDA_API ElementHyperelasticityFEMForceField<sofa::gpu::cuda::CudaVec1fTypes, sofa::geometry::Edge>;
template class ELASTICITY_CUDA_API ElementHyperelasticityFEMForceField<sofa::gpu::cuda::CudaVec2fTypes, sofa::geometry::Edge>;
template class ELASTICITY_CUDA_API ElementHyperelasticityFEMForceField<sofa::gpu::cuda::CudaVec3fTypes, sofa::geometry::Edge>;
template class ELASTICITY_CUDA_API ElementHyperelasticityFEMForceField<sofa::gpu::cuda::CudaVec2fTypes, sofa::geometry::Triangle>;
template class ELASTICITY_CUDA_API ElementHyperelasticityFEMForceField<sofa::gpu::cuda::CudaVec3fTypes, sofa::geometry::Triangle>;
template class ELASTICITY_CUDA_API ElementHyperelasticityFEMForceField<sofa::gpu::cuda::CudaVec2fTypes, sofa::geometry::Quad>;
template class ELASTICITY_CUDA_API ElementHyperelasticityFEMForceField<sofa::gpu::cuda::CudaVec3fTypes, sofa::geometry::Quad>;
template class ELASTICITY_CUDA_API ElementHyperelasticityFEMForceField<sofa::gpu::cuda::CudaVec3fTypes, sofa::geometry::Tetrahedron>;
template class ELASTICITY_CUDA_API ElementHyperelasticityFEMForceField<sofa::gpu::cuda::CudaVec3fTypes, sofa::geometry::Hexahedron>;

namespace cuda
{

void registerCudaElementHyperelasticityFEMForceField(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData("Hyperelasticity on linear beams")
            .add< ElementHyperelasticityFEMForceField<sofa::gpu::cuda::CudaVec1fTypes, sofa::geometry::Edge> >()
            .add< ElementHyperelasticityFEMForceField<sofa::gpu::cuda::CudaVec2fTypes, sofa::geometry::Edge> >()
            .add< ElementHyperelasticityFEMForceField<sofa::gpu::cuda::CudaVec3fTypes, sofa::geometry::Edge> >());

    factory->registerObjects(sofa::core::ObjectRegistrationData("Hyperelasticity on linear triangles")
        .add< ElementHyperelasticityFEMForceField<sofa::gpu::cuda::CudaVec2fTypes, sofa::geometry::Triangle> >()
        .add< ElementHyperelasticityFEMForceField<sofa::gpu::cuda::CudaVec3fTypes, sofa::geometry::Triangle> >());

    factory->registerObjects(sofa::core::ObjectRegistrationData("Hyperelasticity on linear quads")
        .add< ElementHyperelasticityFEMForceField<sofa::gpu::cuda::CudaVec2fTypes, sofa::geometry::Quad> >()
        .add< ElementHyperelasticityFEMForceField<sofa::gpu::cuda::CudaVec3fTypes, sofa::geometry::Quad> >());

    factory->registerObjects(sofa::core::ObjectRegistrationData("Hyperelasticity on linear tetrahedra")
        .add< ElementHyperelasticityFEMForceField<sofa::gpu::cuda::CudaVec3fTypes, sofa::geometry::Tetrahedron> >());

    factory->registerObjects(sofa::core::ObjectRegistrationData("Hyperelasticity on linear hexahedra")
        .add< ElementHyperelasticityFEMForceField<sofa::gpu::cuda::CudaVec3fTypes, sofa::geometry::Hexahedron> >());
}

}


}

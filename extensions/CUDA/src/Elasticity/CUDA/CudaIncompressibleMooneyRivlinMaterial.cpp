#include <Elasticity/component/material/IncompressibleMooneyRivlinMaterial.inl>
#include <Elasticity/CUDA/config.h>
#include <sofa/gpu/cuda/CudaTypes.h>
#include <sofa/core/ObjectFactory.h>

namespace elasticity
{

template class ELASTICITY_CUDA_API IncompressibleMooneyRivlinMaterial<sofa::gpu::cuda::CudaVec1fTypes>;
template class ELASTICITY_CUDA_API IncompressibleMooneyRivlinMaterial<sofa::gpu::cuda::CudaVec2fTypes>;
template class ELASTICITY_CUDA_API IncompressibleMooneyRivlinMaterial<sofa::gpu::cuda::CudaVec3fTypes>;

namespace cuda
{

void registerCudaIncompressibleMooneyRivlinMaterial(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData("Incompressible Mooney-Rivlin material")
        .add< IncompressibleMooneyRivlinMaterial<sofa::gpu::cuda::CudaVec1fTypes> >()
        .add< IncompressibleMooneyRivlinMaterial<sofa::gpu::cuda::CudaVec2fTypes> >()
        .add< IncompressibleMooneyRivlinMaterial<sofa::gpu::cuda::CudaVec3fTypes> >());
}

}

}

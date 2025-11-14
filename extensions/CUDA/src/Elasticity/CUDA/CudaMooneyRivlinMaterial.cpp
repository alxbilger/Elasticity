#include <Elasticity/component/material/MooneyRivlinMaterial.inl>
#include <Elasticity/CUDA/config.h>
#include <sofa/gpu/cuda/CudaTypes.h>
#include <sofa/core/ObjectFactory.h>

namespace elasticity
{

template class ELASTICITY_CUDA_API MooneyRivlinMaterial<sofa::gpu::cuda::CudaVec1fTypes>;
template class ELASTICITY_CUDA_API MooneyRivlinMaterial<sofa::gpu::cuda::CudaVec2fTypes>;
template class ELASTICITY_CUDA_API MooneyRivlinMaterial<sofa::gpu::cuda::CudaVec3fTypes>;

namespace cuda
{

void registerCudaMooneyRivlinMaterial(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData("Mooney-Rivlin material")
        .add< MooneyRivlinMaterial<sofa::gpu::cuda::CudaVec1fTypes> >()
        .add< MooneyRivlinMaterial<sofa::gpu::cuda::CudaVec2fTypes> >()
        .add< MooneyRivlinMaterial<sofa::gpu::cuda::CudaVec3fTypes> >());
}

}

}

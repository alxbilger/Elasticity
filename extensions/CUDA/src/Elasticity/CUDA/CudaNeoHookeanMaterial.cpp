#include <Elasticity/component/material/NeoHookeanMaterial.inl>
#include <Elasticity/CUDA/config.h>
#include <sofa/gpu/cuda/CudaTypes.h>
#include <sofa/core/ObjectFactory.h>

namespace elasticity
{

template class ELASTICITY_CUDA_API NeoHookeanMaterial<sofa::gpu::cuda::CudaVec1fTypes>;
template class ELASTICITY_CUDA_API NeoHookeanMaterial<sofa::gpu::cuda::CudaVec2fTypes>;
template class ELASTICITY_CUDA_API NeoHookeanMaterial<sofa::gpu::cuda::CudaVec3fTypes>;

namespace cuda
{

void registerCudaNeoHookeanMaterial(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData("Neo-Hookean material")
        .add< NeoHookeanMaterial<sofa::gpu::cuda::CudaVec1fTypes> >()
        .add< NeoHookeanMaterial<sofa::gpu::cuda::CudaVec2fTypes> >()
        .add< NeoHookeanMaterial<sofa::gpu::cuda::CudaVec3fTypes> >());
}

}

}

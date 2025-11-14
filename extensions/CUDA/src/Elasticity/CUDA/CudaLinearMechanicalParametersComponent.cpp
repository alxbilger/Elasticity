#include <Elasticity/component/LinearMechanicalParametersComponent.inl>
#include <Elasticity/CUDA/config.h>
#include <sofa/gpu/cuda/CudaTypes.h>

namespace elasticity
{

template class ELASTICITY_CUDA_API LinearMechanicalParametersComponent<sofa::gpu::cuda::CudaVec3fTypes>;
template class ELASTICITY_CUDA_API LinearMechanicalParametersComponent<sofa::gpu::cuda::CudaVec2fTypes>;
template class ELASTICITY_CUDA_API LinearMechanicalParametersComponent<sofa::gpu::cuda::CudaVec1fTypes>;

}

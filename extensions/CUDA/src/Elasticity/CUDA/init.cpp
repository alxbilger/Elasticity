#include <Elasticity/CUDA/init.h>
#include <Elasticity/init.h>
#include <SofaCUDA/init.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/system/PluginManager.h>

namespace elasticity::cuda
{

extern void registerCudaElementCorotationalFEMForceField(sofa::core::ObjectFactory* factory);
extern void registerCudaElementHyperelasticityFEMForceField(sofa::core::ObjectFactory* factory);
extern void registerCudaElementLinearSmallStrainFEMForceField(sofa::core::ObjectFactory* factory);
extern void registerCudaIncompressibleMooneyRivlinMaterial(sofa::core::ObjectFactory* factory);
extern void registerCudaMooneyRivlinMaterial(sofa::core::ObjectFactory* factory);
extern void registerCudaNeoHookeanMaterial(sofa::core::ObjectFactory* factory);
extern void registerCudaStVenantKirchhoffMaterial(sofa::core::ObjectFactory* factory);

void initializePlugin() 
{
    static bool first = true;
    if (first) {
        first = false;

        // make sure that this plugin is registered into the PluginManager
        sofa::helper::system::PluginManager::getInstance().registerPlugin(MODULE_NAME);

        elasticity::initializePlugin();
        sofa::gpu::cuda::init();
    }
}

}

extern "C" 
{
    ELASTICITY_CUDA_API void initExternalModule() 
    {
        elasticity::cuda::initializePlugin();
    }

    ELASTICITY_CUDA_API const char* getModuleName() 
    {
        return elasticity::cuda::MODULE_NAME;
    }

    ELASTICITY_CUDA_API const char* getModuleVersion() 
    {
        return elasticity::cuda::MODULE_VERSION;
    }

    ELASTICITY_CUDA_API const char* getModuleLicense() 
    {
        return "LGPL";
    }

    ELASTICITY_CUDA_API const char* getModuleDescription() 
    {
        return "SOFA plugin for Elasticity.CUDA";
    }

    ELASTICITY_CUDA_API void registerObjects(sofa::core::ObjectFactory* factory)
    {
        elasticity::cuda::registerCudaElementCorotationalFEMForceField(factory);
        elasticity::cuda::registerCudaElementHyperelasticityFEMForceField(factory);
        elasticity::cuda::registerCudaElementLinearSmallStrainFEMForceField(factory);
        elasticity::cuda::registerCudaIncompressibleMooneyRivlinMaterial(factory);
        elasticity::cuda::registerCudaMooneyRivlinMaterial(factory);
        elasticity::cuda::registerCudaNeoHookeanMaterial(factory);
        elasticity::cuda::registerCudaStVenantKirchhoffMaterial(factory);
    }
}

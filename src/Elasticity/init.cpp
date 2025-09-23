#include <Elasticity/init.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/system/PluginManager.h>

namespace elasticity
{

extern void registerElementLinearSmallStrainFEMForceField(sofa::core::ObjectFactory* factory);
extern void registerLinearSmallStrainFEMForceField(sofa::core::ObjectFactory* factory);
extern void registerCorotationalFEMForceField(sofa::core::ObjectFactory* factory);
extern void registerHyperelasticityFEMForceField(sofa::core::ObjectFactory* factory);

void initializePlugin() 
{
    static bool first = true;
    if (first)
    {
        first = false;
        sofa::helper::system::PluginManager::getInstance().registerPlugin(MODULE_NAME);
    }
}

}

extern "C" 
{
    ELASTICITY_API void initExternalModule() 
    {
        elasticity::initializePlugin();
    }

    ELASTICITY_API const char* getModuleName() 
    {
        return elasticity::MODULE_NAME;
    }

    ELASTICITY_API const char* getModuleVersion() 
    {
        return elasticity::MODULE_VERSION;
    }

    ELASTICITY_API const char* getModuleLicense() 
    {
        return "LGPL";
    }

    ELASTICITY_API const char* getModuleDescription() 
    {
        return "SOFA plugin for Elasticity";
    }

    ELASTICITY_API void registerObjects(sofa::core::ObjectFactory* factory)
    {
        elasticity::registerElementLinearSmallStrainFEMForceField(factory);
        elasticity::registerLinearSmallStrainFEMForceField(factory);
        elasticity::registerCorotationalFEMForceField(factory);
        elasticity::registerHyperelasticityFEMForceField(factory);
    }
}

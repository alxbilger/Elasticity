#include <benchmark/benchmark.h>

#include <Elasticity/init.h>

#include <sofa/simulation/common/init.h>
#include <sofa/simulation/graph/init.h>
#include <sofa/simulation/init.h>
#include <sofa/component/init.h>

int main(int argc, char** argv)
{
    char arg0_default[]="benchmark";
    char* args_default=arg0_default;
    if(!argv)
    {
        argc=1;
        argv=&args_default;
    }
    ::benchmark::Initialize(&argc, argv);
    if(::benchmark::ReportUnrecognizedArguments(argc, argv)) return 1;

    sofa::simulation::graph::init();
    sofa::component::init();
    elasticity::initializePlugin();
    sofa::simulation::core::init();

    ::benchmark::RunSpecifiedBenchmarks();
    ::benchmark::Shutdown();

    sofa::simulation::common::cleanup();
    sofa::simulation::graph::cleanup();
    return 0;
}

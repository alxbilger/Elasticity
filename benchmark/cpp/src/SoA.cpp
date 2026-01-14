#include <Elasticity/component/ElementHyperelasticityFEMForceField.h>
#include <Elasticity/impl/MatSoA.h>
#include <Elasticity/impl/VecSoA.h>
#include <benchmark/benchmark.h>

static void BM_Vec3SoA(benchmark::State& state)
{
    elasticity::VecSoA<3, SReal> a(state.range(0));
    elasticity::VecSoA<3, SReal> b(state.range(0));
    elasticity::VecSoA<3, SReal> c(state.range(0));

    for (auto _ : state)
    {
        elasticity::VecSoA<3, SReal>::Add(c, a, b);
    }
}

static void BM_Vec3AoS(benchmark::State& state)
{
    sofa::type::vector<sofa::type::Vec<3, SReal> > a(state.range(0));
    sofa::type::vector<sofa::type::Vec<3, SReal> > b(state.range(0));
    sofa::type::vector<sofa::type::Vec<3, SReal> > c(state.range(0));

    for (auto _ : state)
    {
        for (auto i = 0; i < state.range(0); ++i)
        {
            c[i] = a[i] + b[i];
        }
    }
}

BENCHMARK(BM_Vec3SoA)->RangeMultiplier(10)->Range(1e3, 1e7);
BENCHMARK(BM_Vec3AoS)->RangeMultiplier(10)->Range(1e3, 1e7);

static void BM_Mat3Vec3SoA(benchmark::State& state)
{
    elasticity::VecSoA<3, SReal> a(state.range(0));
    elasticity::MatSoA<3, 3, SReal> b(state.range(0));
    elasticity::VecSoA<3, SReal> c(state.range(0));

    for (auto _ : state)
    {
        elasticity::matrixVectorProduct(a, b, c);
    }
}

static void BM_Mat3Vec3AoS(benchmark::State& state)
{
    sofa::type::vector<sofa::type::Vec<3, SReal> > a(state.range(0));
    sofa::type::vector<sofa::type::Mat<3, 3, SReal> > b(state.range(0));
    sofa::type::vector<sofa::type::Vec<3, SReal> > c(state.range(0));

    for (auto _ : state)
    {
        for (auto i = 0; i < state.range(0); ++i)
        {
            a[i] = b[i] * a[i];
        }
    }
}

BENCHMARK(BM_Mat3Vec3SoA)->RangeMultiplier(10)->Range(1e3, 1e7);
BENCHMARK(BM_Mat3Vec3AoS)->RangeMultiplier(10)->Range(1e3, 1e7);

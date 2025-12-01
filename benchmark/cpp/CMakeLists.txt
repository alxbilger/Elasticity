cmake_minimum_required(VERSION 3.14)
project(ElasticityBenchmark)

include(FetchContent)

# By default, google benchmark requires gtest as a dependency
# The gtest sources in Sofa.Testing are not compatible with google benchmark
# because of the way gtest is found in the cmake script
# Tests are disabled, but they can be enabled by providing either:
# - BENCHMARK_DOWNLOAD_DEPENDENCIES set to true
# - GOOGLETEST_PATH to an existing directory where gtest has been checked out
set(BENCHMARK_ENABLE_GTEST_TESTS FALSE)
FetchContent_Declare(
    googlebenchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG main
)
FetchContent_MakeAvailable(googlebenchmark)

find_package(Sofa.Component REQUIRED)
find_package(Sofa.Simulation.Graph REQUIRED)
sofa_find_package(MultiThreading REQUIRED)
sofa_find_package(Elasticity REQUIRED)

set(SOFABENCHMARK_SRC "src")
set(HEADER_FILES "")
set(SOURCE_FILES ${SOFABENCHMARK_SRC}/main.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES} ${HEADER_FILES})
target_link_libraries(${PROJECT_NAME} PUBLIC benchmark::benchmark)
target_link_libraries(${PROJECT_NAME} PUBLIC Elasticity Sofa.Component Sofa.Simulation.Graph MultiThreading)
target_include_directories(${PROJECT_NAME} PUBLIC ${SOFABENCHMARK_SRC})

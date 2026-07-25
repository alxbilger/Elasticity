# CMakeLists.txt for Elasticity plugin
cmake_minimum_required(VERSION 3.12)
project(Elasticity VERSION 1.0 LANGUAGES CXX)

find_package(Sofa.Core REQUIRED)
sofa_find_package(Sofa.Component.SolidMechanics.FEM.HyperElastic REQUIRED)

set(HEADER_FILES
    src/Elasticity/config.h.in
    src/Elasticity/init.h

    src/Elasticity/component/SoAElementLinearSmallStrainFEMForceField.h
    src/Elasticity/component/SoAElementLinearSmallStrainFEMForceField.inl

    src/Elasticity/impl/IsotropicElasticityTensor.h
    src/Elasticity/impl/MatSoA.h
    src/Elasticity/impl/StrainDisplacement.h
    src/Elasticity/impl/VecSoA.h
    src/Elasticity/impl/VonMisesStressContainer.h
)

set(SOURCE_FILES
    src/Elasticity/init.cpp

    src/Elasticity/component/SoAElementLinearSmallStrainFEMForceField.cpp

)

add_library(${PROJECT_NAME} SHARED README.md ${HEADER_FILES} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} PUBLIC Sofa.Core Sofa.Component.SolidMechanics.FEM.HyperElastic)

sofa_create_package_with_targets(
    PACKAGE_NAME ${PROJECT_NAME}
    PACKAGE_VERSION ${PROJECT_VERSION}
    TARGETS ${PROJECT_NAME} AUTO_SET_TARGET_PROPERTIES
    INCLUDE_SOURCE_DIR "src"
    INCLUDE_INSTALL_DIR ${PROJECT_NAME}
    RELOCATABLE "plugins"
    )
    
cmake_dependent_option(ELASTICITY_BUILD_TESTS "Compile the automatic tests" ON "SOFA_BUILD_TESTS OR NOT DEFINED SOFA_BUILD_TESTS" OFF)
if(ELASTICITY_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

sofa_add_subdirectory(plugin extensions/CUDA Elasticity.CUDA)
sofa_add_subdirectory(application benchmark/cpp ElasticityBenchmark)

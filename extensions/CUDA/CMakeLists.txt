# CMakeLists.txt for Elasticity.CUDA plugin
cmake_minimum_required(VERSION 3.12)
project(Elasticity.CUDA VERSION 1.0 LANGUAGES CXX CUDA)

find_package(Elasticity REQUIRED)
sofa_find_package(SofaCUDA REQUIRED)

set(HEADER_FILES
    src/Elasticity/CUDA/config.h.in
    src/Elasticity/CUDA/init.h
)

set(SOURCE_FILES
    src/Elasticity/CUDA/init.cpp

    src/Elasticity/CUDA/CudaElementCorotationalFEMForceField.cpp
    src/Elasticity/CUDA/CudaElementHyperelasticityFEMForceField.cpp
    src/Elasticity/CUDA/CudaElementLinearSmallStrainFEMForceField.cpp
    src/Elasticity/CUDA/CudaIncompressibleMooneyRivlinMaterial.cpp
    src/Elasticity/CUDA/CudaLinearMechanicalParametersComponent.cpp
    src/Elasticity/CUDA/CudaMooneyRivlinMaterial.cpp
    src/Elasticity/CUDA/CudaNeoHookeanMaterial.cpp
    src/Elasticity/CUDA/CudaStVenantKirchhoffMaterial.cpp
)

add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} PUBLIC Elasticity SofaCUDA)

sofa_create_package_with_targets(
    PACKAGE_NAME ${PROJECT_NAME}
    PACKAGE_VERSION ${PROJECT_VERSION}
    TARGETS ${PROJECT_NAME} AUTO_SET_TARGET_PROPERTIES
    INCLUDE_SOURCE_DIR "src"
    INCLUDE_INSTALL_DIR ${PROJECT_NAME}
    RELOCATABLE "plugins"
    )

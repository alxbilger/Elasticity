# CMakeLists.txt for Elasticity plugin
cmake_minimum_required(VERSION 3.12)
project(Elasticity VERSION 1.0 LANGUAGES CXX)

find_package(Sofa.Core REQUIRED)

set(HEADER_FILES
    src/Elasticity/config.h.in
    src/Elasticity/init.h

    src/Elasticity/FiniteElement.h
    src/Elasticity/FiniteElement[Edge].h
    src/Elasticity/FiniteElement[Hexahedron].h
    src/Elasticity/FiniteElement[Quad].h
    src/Elasticity/FiniteElement[Tetrahedron].h
    src/Elasticity/FiniteElement[Triangle].h

    src/Elasticity/CorotationalFEM.h
    src/Elasticity/CorotationalFEM.inl
    src/Elasticity/CorotationalFEMForceField.h
    src/Elasticity/CorotationalFEMForceField.inl
    src/Elasticity/ElementLinearSmallStrainFEMForceField.h
    src/Elasticity/ElementLinearSmallStrainFEMForceField.inl
    src/Elasticity/LinearFEM.h
    src/Elasticity/LinearFEM.inl
    src/Elasticity/LinearSmallStrainFEMForceField.h
    src/Elasticity/LinearSmallStrainFEMForceField.inl
    src/Elasticity/MatrixTools.h
    src/Elasticity/VectorTools.h
)

set(SOURCE_FILES
    src/Elasticity/init.cpp

    src/Elasticity/CorotationalFEM.cpp
    src/Elasticity/CorotationalFEMForceField.cpp
    src/Elasticity/ElementLinearSmallStrainFEMForceField.cpp
    src/Elasticity/LinearFEM.cpp
    src/Elasticity/LinearSmallStrainFEMForceField.cpp
)

add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} PUBLIC Sofa.Core)

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

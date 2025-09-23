# CMakeLists.txt for Elasticity plugin
cmake_minimum_required(VERSION 3.12)
project(Elasticity VERSION 1.0 LANGUAGES CXX)

find_package(Sofa.Core REQUIRED)

set(HEADER_FILES
    src/Elasticity/config.h.in
    src/Elasticity/init.h

    src/Elasticity/finiteelement/FiniteElement.h
    src/Elasticity/finiteelement/FiniteElement[all].h
    src/Elasticity/finiteelement/FiniteElement[Edge].h
    src/Elasticity/finiteelement/FiniteElement[Hexahedron].h
    src/Elasticity/finiteelement/FiniteElement[Quad].h
    src/Elasticity/finiteelement/FiniteElement[Tetrahedron].h
    src/Elasticity/finiteelement/FiniteElement[Triangle].h

    src/Elasticity/component/BaseElasticityFEMForceField.h
    src/Elasticity/component/BaseElasticityFEMForceField.inl
    src/Elasticity/component/CorotationalFEMForceField.h
    src/Elasticity/component/CorotationalFEMForceField.inl
    src/Elasticity/component/ElementLinearSmallStrainFEMForceField.h
    src/Elasticity/component/ElementLinearSmallStrainFEMForceField.inl
    src/Elasticity/component/HyperelasticityFEMForceField.h
    src/Elasticity/component/HyperelasticityFEMForceField.inl
    src/Elasticity/component/LinearSmallStrainFEMForceField.h
    src/Elasticity/component/LinearSmallStrainFEMForceField.inl

    src/Elasticity/impl/CorotationalFEM.h
    src/Elasticity/impl/CorotationalFEM.inl
    src/Elasticity/impl/ElasticityTensor.h
    src/Elasticity/impl/ElementStiffnessMatrix.h
    src/Elasticity/impl/LinearFEM[all].h
    src/Elasticity/impl/LinearFEM.h
    src/Elasticity/impl/LinearFEM.inl
    src/Elasticity/impl/MatrixTools.h
    src/Elasticity/impl/StrainDisplacement.h
    src/Elasticity/impl/SymmetricTensor.h
    src/Elasticity/impl/VectorTools.h
    src/Elasticity/impl/VecView.h
    src/Elasticity/impl/VonMisesStressContainer.h
)

set(SOURCE_FILES
    src/Elasticity/init.cpp

    src/Elasticity/component/BaseElasticityFEMForceField.cpp
    src/Elasticity/component/CorotationalFEMForceField.cpp
    src/Elasticity/component/ElementLinearSmallStrainFEMForceField.cpp
    src/Elasticity/component/HyperelasticityFEMForceField.cpp
    src/Elasticity/component/LinearSmallStrainFEMForceField.cpp

    src/Elasticity/impl/CorotationalFEM.cpp
    src/Elasticity/impl/LinearFEM.cpp
)

add_library(${PROJECT_NAME} SHARED README.md ${HEADER_FILES} ${SOURCE_FILES})
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

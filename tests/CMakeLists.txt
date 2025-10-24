cmake_minimum_required(VERSION 3.12)
project(Elasticity_test VERSION 1.0)
find_package(Sofa.Testing REQUIRED)

set(SOURCE_FILES
    CorotationalFEM_test.cpp
    FiniteElement_test.cpp
    TetrahedronLinearSmallStrainFEMForceField_test.cpp
    HexahedronLinearSmallStrainFEMForceField_test.cpp
    Material_test.cpp
    MatrixTools_test.cpp
)

sofa_find_package(Sofa.Component.SolidMechanics.Testing REQUIRED)
sofa_find_package(Sofa.Component.Topology.Container.Constant REQUIRED)
add_executable(${PROJECT_NAME} ${SOURCE_FILES})

target_link_libraries(${PROJECT_NAME} Sofa.Testing Elasticity
    Sofa.Component.SolidMechanics.Testing
    Sofa.Component.Topology.Container.Constant)

add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})

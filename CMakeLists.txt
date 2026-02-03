# CMakeLists.txt for Elasticity plugin
cmake_minimum_required(VERSION 3.12)
project(Elasticity VERSION 1.0 LANGUAGES CXX)

find_package(Sofa.Core REQUIRED)
sofa_find_package(Sofa.Component.SolidMechanics.FEM.Elastic REQUIRED)

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

    src/Elasticity/component/BaseElementLinearFEMForceField.h
    src/Elasticity/component/BaseElementLinearFEMForceField.inl
    src/Elasticity/component/ElementCorotationalFEMForceField.h
    src/Elasticity/component/ElementCorotationalFEMForceField.inl
    src/Elasticity/component/ElementHyperelasticityFEMForceField.h
    src/Elasticity/component/ElementHyperelasticityFEMForceField.inl
    src/Elasticity/component/ElementLinearSmallStrainFEMForceField.h
    src/Elasticity/component/ElementLinearSmallStrainFEMForceField.inl
    src/Elasticity/component/ElementMass.h
    src/Elasticity/component/ElementMass.inl
    src/Elasticity/component/ElementPrefab.h
    src/Elasticity/component/ElementPrefab.inl
    src/Elasticity/component/FEMForceField.h
    src/Elasticity/component/FEMForceField.inl
    src/Elasticity/component/HyperelasticMaterial.h
    src/Elasticity/component/HyperelasticMaterial.inl
    src/Elasticity/component/LinearMechanicalParametersComponent.h
    src/Elasticity/component/LinearMechanicalParametersComponent.inl
    src/Elasticity/component/PK2HyperelasticMaterial.h
    src/Elasticity/component/PK2HyperelasticMaterial.inl
    src/Elasticity/component/SoAElementLinearSmallStrainFEMForceField.h
    src/Elasticity/component/SoAElementLinearSmallStrainFEMForceField.inl

    src/Elasticity/component/material/IncompressibleMooneyRivlinMaterial.h
    src/Elasticity/component/material/IncompressibleMooneyRivlinMaterial.inl
    src/Elasticity/component/material/MooneyRivlinMaterial.h
    src/Elasticity/component/material/MooneyRivlinMaterial.inl
    src/Elasticity/component/material/NeoHookeanMaterial.h
    src/Elasticity/component/material/NeoHookeanMaterial.inl
    src/Elasticity/component/material/OgdenMaterial.h
    src/Elasticity/component/material/OgdenMaterial.inl
    src/Elasticity/component/material/StVenantKirchhoffMaterial.h
    src/Elasticity/component/material/StVenantKirchhoffMaterial.inl

    src/Elasticity/component/prefab/CorotationalFEMForceField.h
    src/Elasticity/component/prefab/CorotationalFEMForceField.inl
    src/Elasticity/component/prefab/HyperelasticityFEMForceField.h
    src/Elasticity/component/prefab/HyperelasticityFEMForceField.inl
    src/Elasticity/component/prefab/LinearSmallStrainFEMForceField.h
    src/Elasticity/component/prefab/LinearSmallStrainFEMForceField.inl

    src/Elasticity/impl/ComputeStrategy.h
    src/Elasticity/impl/ElementStiffnessMatrix.h
    src/Elasticity/impl/FullySymmetric4Tensor.h
    src/Elasticity/impl/IdentityMatrix.h
    src/Elasticity/impl/IsotropicElasticityTensor.h
    src/Elasticity/impl/KroneckerDelta.h
    src/Elasticity/impl/LameParameters.h
    src/Elasticity/impl/MajorSymmetric4Tensor.h
    src/Elasticity/impl/MatrixTools.h
    src/Elasticity/impl/MatSoA.h
    src/Elasticity/impl/Strain.h
    src/Elasticity/impl/StrainDisplacement.h
    src/Elasticity/impl/SymmetricTensor.h
    src/Elasticity/impl/VecSoA.h
    src/Elasticity/impl/VecView.h
    src/Elasticity/impl/VectorTools.h
    src/Elasticity/impl/VoigtNotation.h
    src/Elasticity/impl/VonMisesStressContainer.h
    src/Elasticity/impl/trait.h

    src/Elasticity/impl/rotations/HexahedronRotation.h
    src/Elasticity/impl/rotations/IdentityRotation.h
    src/Elasticity/impl/rotations/PolarDecomposition.h
    src/Elasticity/impl/rotations/RotationMethodsContainer.h
    src/Elasticity/impl/rotations/StablePolarDecomposition.h
    src/Elasticity/impl/rotations/TriangleRotation.h
)

set(SOURCE_FILES
    src/Elasticity/init.cpp

    src/Elasticity/component/BaseElementLinearFEMForceField.cpp
    src/Elasticity/component/ElementCorotationalFEMForceField.cpp
    src/Elasticity/component/ElementHyperelasticityFEMForceField.inl
    src/Elasticity/component/ElementLinearSmallStrainFEMForceField.cpp
    src/Elasticity/component/ElementMass.cpp
    src/Elasticity/component/ElementHyperelasticityFEMForceField.cpp
    src/Elasticity/component/ElementPrefab.cpp
    src/Elasticity/component/FEMForceField.cpp
    src/Elasticity/component/HyperelasticMaterial.cpp
    src/Elasticity/component/LinearMechanicalParametersComponent.cpp
    src/Elasticity/component/PK2HyperelasticMaterial.cpp
    src/Elasticity/component/SoAElementLinearSmallStrainFEMForceField.cpp

    src/Elasticity/component/material/IncompressibleMooneyRivlinMaterial.cpp
    src/Elasticity/component/material/MooneyRivlinMaterial.cpp
    src/Elasticity/component/material/OgdenMaterial.cpp
    src/Elasticity/component/material/NeoHookeanMaterial.cpp
    src/Elasticity/component/material/StVenantKirchhoffMaterial.cpp

    src/Elasticity/component/prefab/CorotationalFEMForceField.cpp
    src/Elasticity/component/prefab/HyperelasticityFEMForceField.cpp
    src/Elasticity/component/prefab/LinearSmallStrainFEMForceField.cpp
)

add_library(${PROJECT_NAME} SHARED README.md ${HEADER_FILES} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} PUBLIC Sofa.Core Sofa.Component.SolidMechanics.FEM.Elastic)

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

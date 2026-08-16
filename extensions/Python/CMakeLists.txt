project(Bindings_Elasticity)

set(SOURCE_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Module_Elasticity.cpp

    ${CMAKE_CURRENT_SOURCE_DIR}/src/Binding_HyperelasticMaterial.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Binding_HyperelasticMaterial.cpp
)

if (NOT TARGET SofaPython3::Plugin)
    find_package(SofaPython3 REQUIRED COMPONENTS Elasticity Bindings.Sofa)
endif()

SP3_add_python_module(
    TARGET       ${PROJECT_NAME}
    PACKAGE      Elasticity
    MODULE       Elasticity
    DESTINATION  .
    SOURCES      ${SOURCE_FILES}
    DEPENDS      SofaPython3::Plugin SofaPython3::Bindings.Sofa Elasticity

)
message("-- SofaPython3 bindings for Elasticity will be created.")

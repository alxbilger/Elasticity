# Elasticity

## Description
A SOFA plugin for the development of a clean FEM formulation.

## Compilation

The plugin does not depend on any other plugin or external library.
The procedure to compile the plugin is the same as for any other plugin and is described in the [SOFA documentation](https://sofa-framework.github.io/doc/plugins/build-a-plugin-from-sources/).

## Possible Roadmap

- Linear elasticity
  - [x] A unique linear formulation for any type of element in any dimension
  - [x] A unique corotational formulation for any type of element in any dimension
  - [x] Computation of von Mises stress for visualization (linear)
  - [x] Computation of von Mises stress for visualization (corotational)
  - [ ] Multiple methods to compute the rotation
  - [ ] Rigorous heterogeneous material
- Hyperlaticity
  - [ ] A unique nonlinear formulation for any type of element in any dimension, with any constitutive equation
  - [ ] Common nonlinear constitutive equation
- Other Features
  - [ ] Support for topological changes
  - [ ] Support for computation of the potential energy
- Performances
  - [x] Benchmark: compare to the existing SOFA force fields
  - [ ] Micro-benchmarking to improve performances
  - [ ] Vectorization
  - [ ] CPU parallelization
  - [ ] GPU parallelization
- Support for mixed mesh
  - [ ] Support for connective elements in this plugin
  - [ ] Support for connective elements in SOFA topologies
- Support for higher-order elements
  - [ ] Support for higher-order elements in this plugin
  - [ ] Support for higher-order elements in SOFA topologies
- V&V
  - [ ] Unit testing
  - [ ] Mechanical tests

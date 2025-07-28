# Elasticity

## Description
A SOFA plugin for the development of a clean FEM formulation.

## Possible Roadmap

- Linear elasticity
  - [x] A unique linear formulation for any type of element in any dimension
  - [x] A unique corotational formulation for any type of element in any dimension
  - [x] Computation of von Mises stress for visualization (linear)
  - [x] Computation of von Mises stress for visualization (corotational)
  - [ ] Rigorous heterogeneous material
- Hyperlaticity
  - [ ] A unique nonlinear formulation for any type of element in any dimension, with any constitutive equation
  - [ ] Common nonlinear constitutive equation
- Performances
  - [ ] Benchmark: compare to the existing SOFA force fields
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

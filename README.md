# Elasticity

## Description
A SOFA plugin for the development of a clean FEM formulation.

Elasticity is a plugin implementing the elastic force fields in SOFA in a generic and modular manner.
In particular, a single piece of code allows generating a component for any type of element, in any dimension. 
On top of that, a generic and element-agnostic component (prefab-like) instantiates the appropriate components based on the elements in the topology. 
The emphasis was placed on the readability of the code, allowing to easily compare the code to the equations.

The plugin supports the following formulations:
- Linear elasticity assuming a small strain
- Linear elasticity using a corotational approach
- Hyperelasticity with different constitutive equations

## Compilation

The plugin does not depend on any other plugin or external library.
The procedure to compile the plugin is the same as for any other plugin and is described in the [SOFA documentation](https://sofa-framework.github.io/doc/plugins/build-a-plugin-from-sources/).

## Hyperelastic materials

The plugin supports the following constitutive equations:

### Saint Venant-Kirchhoff

The strain energy density function is:

$$
\psi = \frac{1}{2} \lambda \, \mathrm{tr} (E)^2 + \mu \, E : E
$$

The second Piola-Kirchhoff stress tensor is:

$$
S = 2 \mu E + \lambda \, \mathrm{tr} (E) \, I
$$

And the elasticity tensor is defined in index notation:

$$
\mathbb{C}_{ijkl} = \mu (\delta_{ik} \delta_{jl} + \delta_{il} \delta_{jk}) + \frac{1}{2} \lambda \delta_{ij} \delta_{kl}
$$

Example in a XML scene:

```xml
...
<StVenantKirchhoffMaterial youngModulus="10000" poissonRatio="0.45"/>
<HyperelasticityFEMForceField name="FEM" topology="@Tetra_topo"/>
...
```


## Possible Roadmap

- Linear elasticity
  - [x] A unique linear formulation for any type of element in any dimension
  - [x] A unique corotational formulation for any type of element in any dimension
  - [x] Computation of von Mises stress for visualization (linear)
  - [x] Computation of von Mises stress for visualization (corotational)
  - [ ] Multiple methods to compute the rotation
  - [ ] Rigorous heterogeneous material
- Hyperlaticity
  - [x] A unique nonlinear formulation for any type of element in any dimension, with any constitutive equation
  - [ ] Common nonlinear constitutive equation
    - [x] Saint Venant-Kirchhoff
    - [ ] Neo-Hookean
    - [ ] Stable Neo-Hookean
    - [ ] Arruda Boyce
    - [ ] Mooney-Rivlin
    - [ ] Ogden
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
- Python bindings
  - [ ] Binding for the derivation of a hyperelastic material
  - [ ] Binding for the automatic differentiation of the stress tensor and tangent stiffness matrix

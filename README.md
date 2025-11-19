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

### Neo-Hookean

The strain energy density function is:

$$
\psi = \frac{\mu}{2} (\mathrm{tr}(C) - d) - \mu \log J + \frac{\lambda}{2} (\log J)^2
$$

The second Piola-Kirchhoff stress tensor is:

$$
S = \mu (I - C^{-1}) + \lambda (\log J) C^{-1}
$$

And the elasticity tensor is defined in index notation:

$$
\mathbb{C}_{ijkl} = (\mu - \lambda \log J)(C^{-1}_{i k} C^{-1}_{l j} + C^{-1}_{i l} C^{-1}_{k j}) + \lambda C^{-1}_{l k} C^{-1}_{i j}
$$

### Nearly-incompressible Mooney-Rivlin

The strain energy density function is:

$$
\psi = \mu_{10} (J^{-2/d} I_C - d) + \mu_{01} (J^{-4/d} II_C - d) + \frac{\kappa}{2} (\log J)^2
$$

The second Piola-Kirchhoff stress tensor is:

$$
S = 2 \mu_{10} J^{-2/dim}(I - \frac{1}{d} I_C C^{-1}) + 2 \mu_{01} J^{-4/dim}(I_C I - C - \frac{2}{d} II_C C^{-1}) + \kappa (\log J) C^{-1}
$$

### Ogden

The strain energy density function is:

$$
\psi = F(J) \frac{\mu}{\alpha^2} \sum_{n=1}^3 \left[\lambda_n^{ \frac{\alpha}{2} }\right] - 3 \frac{\mu}{\alpha^2} + \frac{k_0}{2} \ln^2 (J)
$$

The second Piola-Kirchhoff stress tensor is:

$$
S_{ij} = \frac{\mu}{\alpha} F(J) \left(C^{\frac{\alpha}{2} - 1}\right)_{ij} - \frac{\mu}{3 \alpha} F(J) tr\left(C^{\frac{\alpha}{2}}\right) \left( C^{-1}\right)_{ji} + k_0 ln(J) \left( C^{-1}\right)_{ji}
$$

The elasticity tensor is defined in index notation:

```math
\begin{aligned}
\mathbb{C}_{ijkl} &= 2 \frac{\mu}{\alpha} F(J) \sum_{n=1}^3 \left[\left(\frac{\alpha}{2} - 1\right) \lambda_n^{\frac{\alpha}{2} - 2} V_{n,i} V_{n,j} V_{n,k} V_{n,l} \right] \\
&+ \frac{\mu}{\alpha} F(J) \sum_{n=1}^3 \sum_{m \neq n} \left[\frac{\lambda^{\frac{\alpha}{2} - 1}_{n} - \lambda^{\frac{\alpha}{2} - 1}_{m}}{\lambda_{n} - \lambda_{m}} \left( V_{n,i} V_{m,j} V_{m,k} V_{n,l} + V_{n,i} V_{m,j} V_{n,k} V_{m,l} \right) \right] \\ 
&+ \frac{\mu}{9} F(J) tr\left(C^{\frac{\alpha}{2}}\right) \left(C^{-1}\right)_{ji} \left(C^{-1}\right)_{lk} - \frac{\mu}{3} F(J) \left[\left( C^{\frac{\alpha}{2} - 1} \right)_{ji} \left(C^{-1}\right)_{lk} + \left( C^{\frac{\alpha}{2} - 1} \right)_{kl} \left(C^{-1}\right)_{ji} \right] \\
&+ \frac{\mu}{3 \alpha} F(J) tr \left( C^{\frac{\alpha}{2}} \right) \left[\left(C^{-1}\right)_{jk} \left(C^{-1}\right)_{li} + \left(C^{-1}\right)_{jl} \left(C^{-1}\right)_{ki} \right]
\end{aligned}
```

## CUDA

The plugin has an extension called `Elaticity.CUDA`. It is another SOFA plugin. For the moment, it only adds compatibility with CUDA types. No computation is performed on the GPU.

## Possible Roadmap

- Integration in SOFA
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
    - [x] Neo-Hookean
    - [ ] Stable Neo-Hookean
    - [ ] Arruda Boyce
    - [x] Mooney-Rivlin
    - [x] Ogden
  - [ ] Viscosity
  - [ ] Plasticity
  - [ ] Anisotropy
- Other Features
  - [ ] Support for topological changes
  - [ ] Support for computation of the potential energy
  - [ ] Mass
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
- Python bindings
  - [ ] Binding for the derivation of a hyperelastic material
  - [ ] Binding for the automatic differentiation of the stress tensor and tangent stiffness matrix

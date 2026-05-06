# Verification of 2D Linear Elasticity Using the Method of Manufactured Solutions (MMS)

## Mathematical Problem

### 2D linear elasticity (small strain)

$$
\begin{cases}
\displaystyle \frac{\partial \sigma_{xx}}{\partial x} + \frac{\partial \sigma_{xy}}{\partial y} + f_x(x,y) = 0 
\\
\\
\displaystyle \frac{\partial \sigma_{xy}}{\partial x} + \frac{\partial \sigma_{yy}}{\partial y} + f_y(x,y) = 0
\end{cases}
\quad \forall (x,y) \in [0, L] \times [0, L]
$$

### Constitutive law (plane stress, $\nu = 0$)

For $\nu = 0$, the shear modulus reduces to $G = E/2$, so that $2G = E$.

$$
\sigma_{xx} = E \cdot \varepsilon_{xx} \qquad
\sigma_{yy} = E \cdot \varepsilon_{yy} \qquad
\sigma_{xy} = 2G \cdot \varepsilon_{xy} = E \cdot \varepsilon_{xy}
$$

### Strain-displacement relations

$$
\varepsilon_{xx} = \frac{\partial u_x}{\partial x}, \quad
\varepsilon_{yy} = \frac{\partial u_y}{\partial y}, \quad
\varepsilon_{xy} = \frac{1}{2}\left( \frac{\partial u_x}{\partial y} + \frac{\partial u_y}{\partial x} \right)
$$

### Variables and parameters

- $u_x(x,y)$, $u_y(x,y)$: displacement components
- $x$, $y$: spatial coordinates
- $L$: domain side length
- $E$: Young's modulus
- $f_x(x,y)$, $f_y(x,y)$: body forces per unit area

### Boundary conditions

- Dirichlet conditions on $x = 0$ and $x = L$: $u_x = 0$, $u_y = 0$
- Neumann conditions on $y = 0$ and $y = L$: prescribed traction $\mathbf{t} = \boldsymbol{\sigma} \cdot \mathbf{n}$

---

## Discretization

- Standard Galerkin finite element formulation
- Uniform structured mesh on $[0, L] \times [0, L]$, with node spacing $\Delta x = L/(n_x-1)$, $\Delta y = L/(n_y-1)$
- Two element types are considered: **bilinear quadrilateral elements (Q1)** and **linear triangular elements (P1)**

### Gauss Quadrature Rules

#### 1D Gauss–Legendre rules (used for edge integrals)

Integrals over a reference interval $[-1, 1]$ are approximated by:

$$
\int_{-1}^{1} g(\xi)\,d\xi \approx \sum_{k=1}^{n_g} w_k\, g(\xi_k)
$$

| Rule | $n_g$ | $\xi_k$ | $w_k$ |
|------|--------|----------|--------|
| 1-point (midpoint) | 1 | $0$ | $2$ |
| 2-point | 2 | $\pm \dfrac{1}{\sqrt{3}}$ | $1$ |

#### 2D Gauss–Legendre rule on the reference square $[-1,1]^2$ (used for Q1 elements)

$$
\int_{-1}^{1}\int_{-1}^{1} g(\xi, \eta)\,d\xi\,d\eta
\approx \sum_{i=1}^{n_g}\sum_{j=1}^{n_g} w_i\, w_j\, g(\xi_i, \eta_j)
$$

with the 2-point rule ($n_g = 2$): $\xi_{1,2} = \pm 1/\sqrt{3}$, $w_{1,2} = 1$.

The physical coordinates and the Jacobian are obtained via the isoparametric mapping:

$$
\mathbf{x}(\xi,\eta) = \sum_{a=1}^{4} N_a(\xi,\eta)\,\mathbf{x}_a, \qquad
J = \frac{\partial \mathbf{x}}{\partial (\xi,\eta)}, \qquad
d A = |\det J|\, d\xi\, d\eta
$$

where the Q1 shape functions on the reference element are:

$$
N_1 = \tfrac{1}{4}(1-\xi)(1-\eta),\quad
N_2 = \tfrac{1}{4}(1+\xi)(1-\eta),\quad
N_3 = \tfrac{1}{4}(1+\xi)(1+\eta),\quad
N_4 = \tfrac{1}{4}(1-\xi)(1+\eta)
$$

#### 1-point rule on triangles (used for P1 elements)

For a triangle $T$ with area $A_T$ and centroid $\mathbf{x}_c$:

$$
\int_T g(x,y)\,dA \approx A_T \cdot g(\mathbf{x}_c)
$$

This rule is exact for linear functions (consistent with P1 elements).

---

## Body Force Assembly

The nodal force vector associated with the body force is:

$$
\mathbf{F}_i = \int_\Omega \mathbf{f}(x,y)\, \phi_i(x,y)\, dA
$$

### Q1 quadrilateral elements

Using the 2-point Gauss rule over each quad $Q$ with nodes $\{\mathbf{x}_a\}_{a=1}^4$:

$$
\mathbf{F}_i \mathrel{+}= \sum_{k,l} w_k w_l\, N_i(\xi_k, \eta_l)\, \mathbf{f}(\mathbf{x}(\xi_k,\eta_l))\, |\det J(\xi_k, \eta_l)|
$$

### P1 triangular elements — trapezoidal nodal quadrature

For a structured mesh, the integral $\int_\Omega \mathbf{f}\,\phi_i\,dA$ reduces to a nodal weight times $\mathbf{f}(x_i, y_i)$:

$$
\mathbf{F}_i = \mathbf{f}(x_i, y_i) \cdot w_i
$$

with weights depending on the node position:

| Node position | Weight $w_i$ |
|---------------|--------------|
| Interior | $\Delta x\,\Delta y$ |
| Edge (not corner) | $\dfrac{\Delta x\,\Delta y}{2}$ |
| Corner | $\dfrac{\Delta x\,\Delta y}{4}$ |

---

## Neumann Boundary Conditions (Traction)

On boundaries with outward normal $\mathbf{n}$, the traction is $\mathbf{t} = \boldsymbol{\sigma} \cdot \mathbf{n}$. The nodal force contribution is:

$$
\mathbf{F}_i = \int_{\partial \Omega_N} \mathbf{t}(x,y)\, \phi_i(x,y)\, ds
$$

Nodes at $x = 0$ and $x = L$ lying on Neumann boundaries are subject to Dirichlet conditions and are **excluded** from Neumann assembly.

### Q1 quadrilateral elements — 2-point Gauss rule on each edge

For an edge $[k_0, k_1]$ of length $\ell_e$ on $y = 0$ or $y = L$, with the 1D parametrisation $t \in [0,1]$:

$$
\mathbf{F}_{k_0} \mathrel{+}= \sum_{k} w_k\, N_0(\xi_k)\, \mathbf{t}(\mathbf{x}(\xi_k))\, \frac{\ell_e}{2}, \qquad
\mathbf{F}_{k_1} \mathrel{+}= \sum_{k} w_k\, N_1(\xi_k)\, \mathbf{t}(\mathbf{x}(\xi_k))\, \frac{\ell_e}{2}
$$

where $N_0(\xi) = \tfrac{1}{2}(1-\xi)$, $N_1(\xi) = \tfrac{1}{2}(1+\xi)$, $\xi_{1,2} = \pm 1/\sqrt{3}$, $w_{1,2} = 1$.

### P1 triangular elements — trapezoidal nodal quadrature on edges

On $y = 0$ or $y = L$ (excluding Dirichlet endpoints at $x = 0$, $x = L$):

$$
\mathbf{F}_i = \mathbf{t}(x_i, y_i) \cdot w_i
$$

with weights:

| Node position on edge | Weight $w_i$ |
|-----------------------|--------------|
| Interior edge node | $\Delta x$ |
| Endpoint ($x=0$ or $x=L$) | excluded (Dirichlet) |

---

## Error Measures

- $\mathbf{u}_h$: finite element solution
- $\mathbf{u}_{ex}$: manufactured solution

### $L^2$ norm

#### Q1 elements — 2-point Gauss quadrature over each quad

$$
\| \mathbf{u}_h - \mathbf{u}_{ex} \|_{L^2}^2 \approx
\sum_Q \sum_{k,l} w_k w_l\,
\left\| \mathbf{u}_h(\mathbf{x}(\xi_k,\eta_l)) - \mathbf{u}_{ex}(\mathbf{x}(\xi_k,\eta_l)) \right\|^2
|\det J(\xi_k,\eta_l)|
$$

#### P1 elements — centroid quadrature over each triangle

$$
\| \mathbf{u}_h - \mathbf{u}_{ex} \|_{L^2}^2 \approx
\sum_T A_T\, \left\| \mathbf{u}_h(\mathbf{x}_c) - \mathbf{u}_{ex}(\mathbf{x}_c) \right\|^2
$$

where $\mathbf{u}_h(\mathbf{x}_c)$ is obtained by averaging the three nodal values.

### $H^1$ semi-norm

#### Q1 elements — 2-point Gauss quadrature over each quad

$$
| \mathbf{u}_h - \mathbf{u}_{ex} |_{H^1}^2 \approx
\sum_Q \sum_{k,l} w_k w_l \sum_{\alpha,\beta}
\left( \frac{\partial (u_h^\alpha - u_{ex}^\alpha)}{\partial x_\beta}(\mathbf{x}(\xi_k,\eta_l)) \right)^2
|\det J(\xi_k,\eta_l)|
$$

where gradients are computed via the chain rule $\nabla_\mathbf{x} = J^{-T} \nabla_{(\xi,\eta)}$.

#### P1 elements — element-constant gradient

Since P1 gradients are constant on each triangle, the exact Gauss rule (any interior point) gives:

$$
| \mathbf{u}_h - \mathbf{u}_{ex} |_{H^1}^2 \approx
\sum_T A_T \sum_{\alpha,\beta}
\left( \frac{\partial (u_h^\alpha - u_{ex}^\alpha)}{\partial x_\beta}\bigg|_T \right)^2
$$

evaluated at the triangle centroid for the exact part.

---

## Manufactured Solution

### Displacement field

$$
u_{x,\text{ex}}(x,y) = \frac{x^2 (L - x)}{L^2}
$$

$$
u_{y,\text{ex}}(x,y) = \frac{x (L - x) y}{L^2}
$$

### Derivatives

$$
\frac{\partial u_x}{\partial x} = \frac{2xL - 3x^2}{L^2}, \qquad
\frac{\partial u_x}{\partial y} = 0
$$

$$
\frac{\partial u_y}{\partial x} = \frac{(L - 2x)\,y}{L^2}, \qquad
\frac{\partial u_y}{\partial y} = \frac{x(L-x)}{L^2}
$$

### Stress components ($\nu = 0$)

$$
\sigma_{xx} = E\,\frac{2xL - 3x^2}{L^2}, \qquad
\sigma_{yy} = E\,\frac{x(L-x)}{L^2}, \qquad
\sigma_{xy} = \frac{E}{2}\,\frac{(L - 2x)\,y}{L^2}
$$

### Source terms (body forces)

From equilibrium $\nabla \cdot \boldsymbol{\sigma} + \mathbf{f} = 0$:

$$
f_x(x,y) = -\frac{\partial \sigma_{xx}}{\partial x} - \frac{\partial \sigma_{xy}}{\partial y}
= -\frac{E(2L - 6x)}{L^2} - \frac{E(L - 2x)}{2L^2}
$$

$$
f_y(x,y) = -\frac{\partial \sigma_{xy}}{\partial x} - \frac{\partial \sigma_{yy}}{\partial y}
= \frac{E\,y}{L^2}
$$

### Boundary conditions

#### Dirichlet (imposed displacement) — $x = 0$ and $x = L$

$$
u_x = 0, \qquad u_y = 0
$$

#### Neumann (imposed traction) — $y = 0$ (outward normal $\mathbf{n} = (0,-1)^T$)

$$
t_x = -\sigma_{xy}(x, 0) = 0, \qquad
t_y = -\sigma_{yy}(x, 0) = -E\,\frac{x(L-x)}{L^2}
$$

#### Neumann (imposed traction) — $y = L$ (outward normal $\mathbf{n} = (0,+1)^T$)

$$
t_x = \sigma_{xy}(x, L) = \frac{E}{2}\,\frac{(L - 2x)\,L}{L^2} = \frac{E(L-2x)}{2L}, \qquad
t_y = \sigma_{yy}(x, L) = E\,\frac{x(L-x)}{L^2}
$$

---

## Source Term and Traction Discretization (Summary)

### Q1 elements — 2-point Gauss rule

All integrals (body forces and Neumann tractions) use the 2-point Gauss–Legendre rule:

- **Body forces**: 2D tensor-product rule ($2\times 2$ points per quad element)
- **Neumann tractions**: 1D 2-point rule per boundary edge

### P1 elements — trapezoidal nodal quadrature

**Body forces** (interior nodes, edge nodes not at corners, corner nodes):

$$
\mathbf{F}_i = \mathbf{f}(x_i, y_i) \cdot \Delta x\,\Delta y, \quad
\mathbf{F}_i = \mathbf{f}(x_i, y_i) \cdot \frac{\Delta x\,\Delta y}{2}, \quad
\mathbf{F}_i = \mathbf{f}(x_i, y_i) \cdot \frac{\Delta x\,\Delta y}{4}
$$

**Neumann tractions** on $y = \text{const}$ edges (interior nodes only, endpoints excluded):

$$
\mathbf{F}_i = \mathbf{t}(x_i, y_i) \cdot \Delta x
$$

---

## Expected Convergence Rates

For a smooth manufactured solution on a uniform mesh, the expected asymptotic rates are:

| Element | $L^2$ error | $H^1$ semi-norm |
|---------|-------------|-----------------|
| Q1 quad | $O(h^2)$ | $O(h^1)$ |
| P1 tri  | $O(h^2)$ | $O(h^1)$ |

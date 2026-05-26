# Verification of 1D Linear Elasticity Using the Method of Manufactured Solutions (MMS)

## Mathematical Problem

### 1D linear elasticity (small strain)

$$
E \frac{d^2 u}{dx^2} + f(x) = 0 \quad \forall x \in [0, L]
$$

### Variables and parameters

- $u(x)$: axial displacement  
- $x$: spatial coordinate  
- $L$: domain length  
- $E$: Young's modulus  
- $f(x)$: body force per unit length  

### Boundary conditions

Dirichlet and/or Neumann boundary conditions at $x = 0$ and $x = L$.

## Non-dimensionalization

With $\hat{x} = x/L$, $\hat{u}(\hat{x}) = u(L\hat{x})$, $\hat{E} = E/L$, $\hat{f}(\hat{x}) = f(L\hat{x})$, the problem becomes:

$$
\hat{E}\frac{d^2 \hat{u}}{d\hat{x}^2} + \hat{f}(\hat{x}) = 0 \quad \forall \hat{x} \in [0,1]
$$

From this point on, all symbols refer to the non-dimensional quantities; the hats are dropped.

## Discretization  

- Standard Galerkin finite element formulation
- Linear finite elements (P1)  
- Uniform mesh on $[0, L]$

Integrals over each element $e = [x_i, x_{i+1}]$ of length $h = x_{i+1} - x_i$ are approximated by:

$$
\int_{x_i}^{x_{i+1}} g(x)\,dx \approx \frac{h}{2} \sum_{k=1}^{n_g} w_k\, g(x_k^{(e)})
$$

where the Gauss points in physical coordinates are:

$$
x_k^{(e)} = \frac{x_i + x_{i+1}}{2} + \frac{h}{2}\,\xi_k
$$

The rule is determined by the choice of $n_g$, $\xi_k$, and $w_k$:

| Rule | $n_g$ | $\xi_k$ | $w_k$ |
|------|--------|----------|--------|
| 1-point (midpoint) | 1 | $0$ | $2$ |
| 2-point | 2 | $\pm \frac{1}{\sqrt{3}}$ | $1$ |

### Body Force

The nodal force vector associated with the body force is defined as:

$$
F_i = \int_0^L f(x)\,\phi_i(x)\,dx
$$

where $\phi_i(x)$ are the finite element shape functions.

### Error Measures

- $u_h$: finite element solution  
- $u_{ex}$: manufactured solution

#### L² norm

$$
\| u_h - u_{ex} \|_{L^2}^2 \approx \sum_e \frac{h}{2} \sum_{k=1}^{n_g} w_k\, \big(u_h(x_k^{(e)}) - u_{ex}(x_k^{(e)})\big)^2
$$

evaluated with the 1-point rule.

#### H¹ semi-norm

$$
| u_h - u_{ex} |_{H^1}^2 \approx \sum_e \frac{h}{2} \sum_{k=1}^{n_g} w_k\, \left( \frac{du_h}{dx}(x_k^{(e)}) - \frac{du_{ex}}{dx}(x_k^{(e)}) \right)^2
$$

evaluated with the 2-point rule.

---

## Case 1: Quadratic Function

### Manufactured Solution

$$
u_{ex}(x) = x(1 - x)
$$

### Source Term

$$
f(x) = 2E
$$

### Boundary Conditions

$$
u(0) = 0
$$

$$
\left.\frac{d u}{d x}\right|_{1} = -1
$$

Neumann force:

$$
F_N = -E
$$


### Source Term Discretization

With a constant $f$, the 1-point quadrature rule results to:

$$
F_0 = \frac{h}{2} f, \quad F_i = h\, f, \quad F_N = \frac{h}{2} f
$$

---

## Case 2: Sinusoidal Function

### Manufactured Solution

$$
u_{ex}(x) = \sin(2\pi x)
$$

### Source Term

$$
f(x) = 4 E \pi^2 \sin(2\pi x)
$$

### Boundary Conditions

$$
u(0) = 0
$$

$$
\left.\frac{du}{dx}\right|_{1} = 2\pi \cos(2\pi) = 2\pi
$$

Neumann force:

$$
F_N = 2\pi E
$$

### Source Term Discretization

With a 2-point quadrature rule:

$$
F_i^{(e)} = \frac{h}{2} \sum_{k=1}^{2} w_k\, f(x_k^{(e)})\, \phi_i(x_k^{(e)})
$$

---

## Case 3: Cubic Function

### Manufactured Solution

$$
u_{ex}(x) = x^2 (1 - x)
$$

### Source Term

$$
f(x) = E (6x - 2)
$$

### Boundary Conditions

$$
u(0) = 0
$$

$$
\left.\frac{du}{dx}\right|_{1} = -1
$$

Neumann force:

$$
F_N = -E
$$

### Source Term Discretization

With a 2-point quadrature rule:

$$
F_0 = \frac{h}{6}\left[2f(x_0) + f(x_1)\right], \quad F_i = h\, f(x_i), \quad F_N = \frac{h}{6}\left[f(x_{N-1}) + 2f(x_N)\right]
$$

$$
F_i^{(e)} = \frac{h}{2} \sum_{k=1}^{2} w_k\, f(x_k^{(e)})\, \phi_i(x_k^{(e)})
$$

---
## Case 4: Exponential Function

### Manufactured Solution


$$u_{ex}(x) = e^x - 1$$

### Source Term

$$
f(x)=−E exf(x) = -E\, e^xf(x)=−Eex
$$

### Boundary Conditions


$$ 
u(0)=0
$$


$$
\left.\frac{du}{dx}\right|_{1} = e
$$

Neumann force:

$$
F_N = E\, e
$$


### Source Term Discretization

With a 2-point quadrature rule:

$$
F_i^{(e)} = \frac{h}{2} \sum_{k=1}^{2} w_k\, f(x_k^{(e)})\, \phi_i(x_k^{(e)})
$$
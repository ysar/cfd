

## Euler equations ##

1D Euler equations -

$$
    \frac{\partial \mathbf{U}}{\partial t} + \frac{\partial \mathbf{F}}{\partial x} = 0
$$

where 

$$
    \mathbf{U} = \left[ \begin{array}{c}
    \rho\\
    \rho u\\
    \rho E
    \end{array} \right] 
$$

$$
    \mathbf{F} = \left[ \begin{array}{c}
    \rho u\\
    \rho u^2 + p\\
    u \left(\rho E + p\right)
    \end{array} \right] 
$$

$$
    p = \left( \gamma - 1\right) \left( \rho E - \frac{\rho u^2}{2}\right)
$$

And $\gamma=1.4$ for air. Also important is the use of the sound speed

$$
    a = \sqrt{\frac{\gamma p}{\rho}}
$$

## MHD equations ##

The ideal MHD equations can be given by - 

The 1D ideal MHD equations can be given as - 

## Approximate fluxes at Riemmann interfaces ##

## Discontinuous Galerkin Method ##
### 1D set of conservation laws ###


We seek to find weak solutions to this pde such that for a function \phi,

$$
    \int_\Omega \phi \left( \frac{\partial \mathbf{U}}{\partial t} + \frac{\partial \mathbf{F}}{\partial x} \right) \,dx = 0
$$

In DG method we assume the solution at element k to take the form

$$
    \mathbf{u}_k = \sum_{j=1}^{p+1} \mathbf{U}_{k,j}(t)\, \phi_j(x)
$$

Where \phi_j(x) are the basis functions we use to represent the solution in element $$k$$ and $$\mathbf{U}_{k,j}$$ are the expansion coefficients for element k for the j th order
polynomial.

Rewriting the weak form of the pde for one element k,
$$
    \int_{\Omega_k} \phi_i \left( \frac{\partial \mathbf{U}}{\partial t} + \frac{\partial \mathbf{F}}{\partial x} \right) \,dx = 0
$$

Which can be written as 

$$
    \int_{\Omega_k} \phi_i  \frac{\partial \mathbf{U}}{\partial t}  \,dx - \int_{\Omega_k} \frac{\partial \phi_i}{\partial x} \mathbf{F} \,dx + \left[ \phi \mathbf{F}\right]_{x=k-1/2}^{x=k+1/2} = 0}
$$

Which can be written as 

$$
    \mathbf{M}_k\frac{\partial \mathbf{U}_k}{\partial t} + \mathbf{R}_k = 0
$$

Which is what we will solve for.


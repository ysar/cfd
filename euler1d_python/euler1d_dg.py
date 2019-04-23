#!/usr/bin/env python3
"""
80chars -> 
--------------------------------------------------------------------------------
Solve the 1D Euler equations using the Discontinuous-Galerkin formulation.
There are many moving parts to a DG code. I will try my best to disentangle all
of them

methods we need (in order of - what comes to mind first)
---------------
    - scipy.integrate.fixed_quad()   [for Gaussian quadrature]
    - calc_mass_matrix
    - calc_residual1_integrand
    - calc_residual2

    - calc_basis_func
    - 
    - advance_in_time
    - reconstruct_solution         [May need this at the end to plot solution]
-------------------------------------------------------------------------------
NOTES:
    - To ravel matrix to vector use <matrix>.ravel()
    - To unravel vector to matrix use <vector>.reshape()

"""
import numpy as np
from euler1d_fv import GAMMA, calc_p, calc_F, calc_Fhat


def convert_to_conservative(u):
    """Given vector of primitive variables, convert to conservative"""
    u[2] = u[2] / (GAMMA - 1) + u[0] * u[1]**2 / 2
    u[1] = u[0] * u[1]
    return u

    
class quadrature(object):
    def __init__(self, degree):
        """Class containing quadrature related constants"""
        if degree == 3:
            self.nquad = 3
            self.epsilon = np.array([-np.sqrt(3/5), 0, np.sqrt(3/5)])
            self.weights = np.array([5/9, 8/9, 5/9])

        elif degree == 2:
            self.nquad = 2
            self.epsilon = np.array([-np.sqrt(1/3), np.sqrt(1/3)])
            self.weights = np.array([1.0, 1.0])

        elif degree == 1:
            self.nquad = 1
            self.epsilon = np.array([0.0])
            self.weights = np.array([1.0])
            
def calc_phi(eps, p, i):
    """
    Given epsilon, return phi_i(epsilon)
    epsilon  - [-1, 1] reference space
    p        - Order of accuracy (# of polynomials = p + 1) 
    i        - i'th basis function

    Assumes a Lagrange basis. 
    """
    phi = 1
    epsilon = np.linspace(-1, 1, p + 1)   # Assumes equidistant points in reference space
    for j in range(p + 1):
        if i != j:
            phi *= (eps - epsilon[j]) / (epsilon[i] - epsilon[j])

    return phi


def calc_dphi_depsilon(eps, p, i):
    """
    Calculate the derivative of the legendre basis functions in reference space
    """
    dphideps = 0
    epsilon = np.linspace(-1, 1, p + 1)

    denom = 1
    for j in range(p + 1):
        if i != j:
            denom *= (epsilon[i] - epsilon[j])

    dphideps = 0
    for j in range(p + 1):
        if i != j:
            product = 1
            for k in range(p + 1):
                if k != j:
                    product *= (eps - epsilon[k])
            dphideps += product
            
    return dphideps / denom


def calc_M(p):
    """
    Calculate the mass matrix for a reference element.
    Size = (p+1, p+1)
    If we use phi with degree p then the mass matrix has a degree of 
    2p. So we need to use atleast (p+1) points for the quadrture which can 
    exactly describe a polynomial of degree (2*(p+1) - 1) = (2p + 1).  
    """
    M = np.zeros((p+1, p+1))
    quad = quadrature(p+1)

    # Does not account for the change in reference basis. For unequal
    # grid spacing / 2D problems we need to construct a different mass
    # matrix for each element. 
    for i in range(p+1):
        for j in range(p+1):
            for q in range(quad.nquad):
                M[i, j] += calc_phi(quad.epsilon[q], p, i) * \
                           calc_phi(quad.epsilon[q], p, j) * \
                           quad.weights[q]
    return M

        
def construct_initial_state(N, p, typeinitial='freestream'):
    """
    Algorithm: (if using Lagrange basis)
        1. Construct solution u for N * (p+1) points. 
        2. For each point k
            2.1. U(k, i) = u(k, j) if   i == j    else   0
    """
    # Specify simple initial conditions based on cell averages. (will need to modify in the future)
    uc = np.zeros((N, 3))
    
    if typeinitial == 'freestream':
        uc[:, 0] = 1.0
        uc[:, 1] = 1.0
        uc[:, 2] = 1.0

    elif typeinitial == 'sod':
        midway = int(N/2 if N%2 == 0 else (N+1) / 2)
        uc[:midway, 0] = 1.0
        uc[:midway, 1] = 0.0
        uc[:midway, 2] = 1.0
        uc[midway:, 0] = 0.125
        uc[midway:, 1] = 0.0
        uc[midway:, 2] = 0.1
    
    uc[:, 2] = uc[:, 2] / (GAMMA - 1) + uc[:, 0] * uc[:, 1]**2 / 2
    uc[:, 1] = uc[:, 0] * uc[:, 1]

    # The good thing about a Lagrange basis is that the
    # solution at each Lagrange node depends only on one
    # basis function.
    # i.e. u_{x=node1} = U_1 phi_1(x=node1)
    # but phi_1 (x=node1) == 1
    # therefore ----->    u_{x=node1} = U_1
    # Calculate the expansion coefficients for each cell.
    U = np.zeros((N, p + 1, 3))
    for i in range(p + 1):
        U[:, i, :] = uc[:, :]

    return U


def reconstruct_plot_solution(U, xfaces):
    """
    Reconstruct and plot the solution.
    u_k = sum over all basis functions (U_i phi_i)
    """
    Nrec = 10   # Number of reconstruction points in each cell.
    N = U.shape[0]
    p = U.shape[1] - 1
    s = U.shape[2]

    # Create array of x at specified points.
    # We will calculate * u * at these points
    xrec = np.zeros((N, Nrec))
    for k in range(N):
            xrec[k, :] = np.linspace(xfaces[k], xfaces[k + 1], Nrec)

    print(xrec.shape)
    
    # Calculate u at each point.
    u = np.zeros((N, Nrec, s))
    for k in range(N):
        for i in range(Nrec):
            eps = (xrec[k, i] - xfaces[k]) / (xfaces[k + 1] - xfaces[k]) * 2 - 1
            for m in range(s):
                for j in range(p + 1):
                    u[k, i, m] += U[k, j, m] * calc_phi(eps, p, j)


    u[:, :, 2] = (GAMMA - 1) * (u[:, :, 2] - u[:, :, 1]**2 / u[:, :, 0] / 2)
    u[:, :, 1] = u[:, :, 1] / u[:, :, 0]
    
    # Plot it.
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(3, 1, figsize=(11, 8.5))
    for m in range(s):
        for k in range(N):
            ax[m].plot(xrec[k, :], u[k, :, m])

    ax[0].set_ylabel(r'$\rho$')
    ax[1].set_ylabel('$u$')
    ax[2].set_ylabel('$p$')
    ax[2].set_xlabel('$x$')
    
    plt.show()
    return xrec, u


def calc_residual_part1(U):
    """Calc integral of dphi/dx * F over the element k"""
    N = U.shape[0]
    p = U.shape[1] - 1
    s = U.shape[2]
    quad = quadrature(p+1)
    Res1 = np.zeros((N, p + 1, s))
    
    for k in range(N):

        # Calculate u and then F at quadrature points.
        # Also calculate dphi_depsilon in the meantime.
        uk = np.zeros((quad.nquad, s))
        Fk = np.zeros((quad.nquad, s))
        dphi_deps = np.zeros((quad.nquad, p + 1))
        
        for i in range(quad.nquad):
            eps = quad.epsilon[i]
            for m in range(s):                
                for j in range(p + 1):
                    dphi_deps[i, j] = calc_dphi_depsilon(eps, p, j)
                    uk[i, m] += U[k, j, m] * calc_phi(eps, p, j)
                    
            Fk[i, :] = calc_F(uk[i, :])

            
        # Do the Gaussian integral for each j and s.
        for j in range(p + 1):
            for m in range(s):
                for i in range(quad.nquad):
                    Res1[k, j, m] += dphi_deps[i, j] * Fk[i, m] * quad.weights[i]
                    
    return Res1


def calc_residual_part2(U):
    """Calculate the second part of the residual-
    (phi F)_rightface - (phi F)_leftface 
    """
    flux = 'roe'
    N = U.shape[0]
    p = U.shape[1] - 1
    s = U.shape[2]
    Res2 = np.zeros((N, p + 1, s))
    # u_sod_left = convert_to_conservative(np.array([1.0, 1.0, 1.0]))
    # u_sod_right = convert_to_conservative(np.array([1.0, 1.0, 1.0]))
    
    u_sod_left = convert_to_conservative(np.array([1.0, 0.0, 1.0]))
    u_sod_right = convert_to_conservative(np.array([0.125, 0.0, 0.1]))
        
    for k in range(N):

        if k == 0:              # Boundary condition 1
            
            Fstar_left = calc_F(u_sod_left)
            
            uL_right = U[k, p, :]
            uR_right = U[k + 1, 0, :]
            Fstar_right = calc_Fhat(uL_right, uR_right, typeflux=flux)
            
        elif k == N - 1:        # Boundary condition 2
            uL_left = U[k - 1, p, :]
            uR_left = U[k, 0, :]
            Fstar_left = calc_Fhat(uL_left, uR_left, typeflux=flux)

            Fstar_right = calc_F(u_sod_right)
            
        else:
            # Calculate uL_left , uR_left, uL_right and uR_right
            uL_left = U[k - 1, p, :] 
            uR_left = U[k, 0, :]
            uL_right = U[k, p, :] 
            uR_right = U[k + 1, 0, :]

            Fstar_left = calc_Fhat(uL_left, uR_left, typeflux=flux)        
            Fstar_right = calc_Fhat(uL_right, uR_right, typeflux=flux)

        for j in range(p + 1):
            Res2[k, j, :] = calc_phi(1.0, p, j) * Fstar_right - calc_phi(-1.0, p, j) * Fstar_left 

    return Res2
 


# ==========================================================================    
if __name__ == '__main__':

    N = 51
    p = 2
    xmin = 0
    xmax = 1
    xfaces = np.linspace(xmin, xmax, N + 1)
    U = construct_initial_state(N, p, typeinitial='sod')
    M = calc_M(p) * (xfaces[1] - xfaces[0]) / 2
    invM = np.linalg.inv(M)
    
    # reconstruct_plot_solution(U, xfaces)
    
    t = 0.0
    dt = 0.0001
    tmax = 0.0001
    
    while t < tmax:

        Res1 = calc_residual_part1(U)
        Res2 = calc_residual_part2(U)
        Res = Res2 - Res1
        
        # print(Res)
        
        for k in range(N):
            for s in range(3):
                U[k, :, s] += dt * np.matmul(-invM, Res[k, :, s])

        t += dt
        
    xrc, urc = reconstruct_plot_solution(U, xfaces)
        
    
    
         
    


    



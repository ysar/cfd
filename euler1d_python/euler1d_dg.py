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
"""
import numpy as np
from euler1d_fv import GAMMA, calc_p, calc_F, calc_Fhat

class quadrature(object):
    def __init__(self, degree):
        """Class containing quadrature related constants"""
        if degree == 3:
            self.nquad = 3
            self.epsilon = np.array([-np.sqrt(3/5), 0, np.sqrt(3/5)])
            self.weights = np.array([5/9, 8/9, 5/9])
            

def calc_phi(eps, p, i):
    """
    Given epsilon, return phi_i(epsilon)
    epsilon  - [-1, 1] reference space
    p        - Order of accuracy (# of polynomials = p + 1) 
    i        - i'th basis function

    Assumes a Lagrange basis. 
    """
    phi = 1
    epsilon = np.linspace(-1, 1, p + 1)
    for j in range(p + 1):
        if i != j:
            phi *= (eps - epsilon[j]) / (epsilon[i] - epsilon[j])

    return phi


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

    for i in range(p+1):
        for j in range(p+1):
            M[i, j] = 0
            for q in range(quad.nquad):
                M[i, j] += calc_phi(quad.epsilon[q], p, i) * \
                           calc_phi(quad.epsilon[q], p, j) * \
                           quad.weights[q]
    return M


def construct_initial_state(N, p):
    """
    Algorithm: (if using Lagrange basis)
        1. Construct solution u for N * (p+1) points. 
        2. For each point k
            2.1. U(k, i) = u(k, j) if   i == j    else   0
    """
    midway = N/2 if (N+1)%2 == 0 else (N+1) / 2
    u = np.zeros((3, N, p + 1))
    u[0, :midway, :] = 1.0
    u[1, :midway, :] = 0.0
    u[2, :midway, :] = 1.0
    u[0, midway:, :] = 0.125
    u[1, midway:, :] = 0.0
    u[2, midway:, :] = 0.1
    u[2, :, :] = u[2, :, :] / (GAMMA - 1) + u[0, :, :] * u[1, :, :]**2 / 2
    u[1, :, :] = u[0, :, :] * u[1, :, :]

    U = np.zeros((3, N, p + 1))
    for k in range(N):
        for i in range(p+1):
            U[:, k, i] = u[:, i]
                
    
if __name__ == '__main__':
    print(calc_M(2))
    

    
         
    


    



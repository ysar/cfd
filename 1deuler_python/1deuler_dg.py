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
-------------------------------------------------------------------------------
"""
from 1deuler_fv import GAMMA, calc_p, calc_F, calc_Fhat

    

    
         
    


    



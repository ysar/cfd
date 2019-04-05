#!/usr/bin/env python3
"""
Code to solve the 1D Euler equations
using a finite-volume method.

Algorithm:
1. Create mesh for N cells.
2. For time t = 0 to tmax
    2.1. Calculate fluxes at each of
         the N+1 interfaces.
         (for different choice of fluxes)

    If time-stepping == explicit
       Calculate U_n+1 
    else
       Create U_n+1 matrix. 
       Solve Ax=b system. 
"""
import numpy as np

# CONSTANTS
GAMMA = 1.4

# --------------------------------------------------------
def set_initial_state(N, sampletest=None):
    """Set initial state based on the test number chosen
    (Rho, RhoU, E)"""
    U = np.zeros((3, N))
    midway = int(N/2) if N % 2 == 0 else int(N/2 + 1)
    
    # A bunch of sample tests taken from the Toro book.
    # Rho, U, P
    if sampletest == 1:
        U[0, :midway] = 1.0
        U[1, :midway] = 0.0
        U[2, :midway] = 1.0
        
        U[0, midway:] = 0.125
        U[1, midway:] = 0.0
        U[2, midway:] = 0.1

    elif sampletest == 2:
        U[0, :midway] = 1.0
        U[1, :midway] = -2.0
        U[2, :midway] = 0.4
        
        U[0, midway:] = 1.0
        U[1, midway:] = 2.0
        U[2, midway:] = 0.4

    elif sampletest == 3:
        U[0, :midway] = 1.0
        U[1, :midway] = 0.0
        U[2, :midway] = 1000.0
        
        U[0, midway:] = 1.0
        U[1, midway:] = 0.0
        U[2, midway:] = 0.01

    elif sampletest == 4:
        U[0, :midway] = 1.0
        U[1, :midway] = 0.0
        U[2, :midway] = 0.01
        
        U[0, midway:] = 1.0
        U[1, midway:] = 0.0
        U[2, midway:] = 100.0

    elif sampletest == 5:
        U[0, :midway] = 5.99924
        U[1, :midway] = 19.5975
        U[2, :midway] = 480.894
        
        U[0, midway:] = 5.99242
        U[1, midway:] = -6.19633
        U[2, midway:] = 46.0950
        

    # Convert to conservative variables
    U[2, :] = U[0, :] * (
        U[1, :]**2 / 2 + U[2, :] / (GAMMA - 1) / U[0, :])
    U[1, :] = U[1, :] *  U[0, :]

    return U

# -----------------------------------------------------------
def calculate_cell_F(U):
    """Calculate the conservative flux vector"""
    F = np.zeros(U.shape)

    p = (U[2, :] / U[1, :] - (
        U[1, :] / U[0, :])**2 / 2) * (GAMMA - 1) * U[0, :]  
    
    F[0, :] = U[1, :]
    F[1, :] = U[1, :]**2 / U[0, :] + p
    F[2, :] = U[1, :] / U[0, :] * (U[2, :] + p)

    return F
    
# -----------------------------------------------------------
def calculate_fluxes(x, xface, U, F, typeflux=None):
    """
    Use the conventional fluxes to calculate the flux
    values at the faces.
    """
    N = len(faces)
    fluxes = np.zeros((3, N + 1))

    if typeflux == 'Rusanov':
        for i in range(N + 1):
            fluxes[:, i] = 
    
# ----------------------------------------------------------
if __name__ == '__main__':

    N = 100

    # Inputs
    xmin = 0
    xmax = 1
    tmax = 1.0
    dt = 0.001
    flux = 'Rusanov'
    
    
    # Discretize domain
    xface = np.linspace(xmin, xmax, N+1)
    x = np.zeros(N)
    for i in range(N):
        x[i] = 0.5 * (faces[i] + faces[i+1])
    dx = np.abs(x[1] - x[0])
    
    # Initialize U at t=0
    U = set_initial_state(sampletest=1)
    
    # Main loop
    t = 0.0
    while t <= tmax:
        F = calculate_cell_F(U)
        fluxes = calculate_fluxes(x, xface, U, F, typeflux=flux)
        U = advance(U, fluxes, deltaT, x)

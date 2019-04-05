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
import matplotlib.pyplot as plt

# CONSTANTS
GAMMA = 1.4


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
    U[2, :] = 0.5 * U[0, :] * U[1, :]**2 + U[2, :] / (GAMMA - 1)
    U[1, :] = U[1, :] *  U[0, :]

    return U


def calculate_cell_F(U):
    """Calculate the conservative flux vector"""
    F = np.zeros(3)

    p = (GAMMA - 1) * (U[2] - 0.5 * U[1]**2 / U[0])
    
    F[0] = U[1]
    F[1] = U[1]**2 / U[0] + p
    F[2] = U[1] / U[0] * (U[2] + p)

    return F


def calculate_fluxes(U, F, N, typeflux=None):
    """
    Use the conventional fluxes to calculate the flux
    values at the faces.

    Ignores the boundary faces. i.e. iface=0 and iface=N+1.
    BCs are specified in another method. 
    """
    fluxes = np.zeros((3, N + 1))
    
    if typeflux == 'Rusanov':

        p = (GAMMA - 1) * (U[2, :] - 0.5 * U[1, :]**2 / U[0, :])
        smax = np.abs(U[1, :] / U[0, :]) + np.sqrt(GAMMA * p / U[0, :])
        
        for i in range(1, N):
            fluxes[:, i] = 0.5 * (F[:, i - 1] + F[:, i]) \
                - 0.5 * smax[i] * (U[:, i] - U[:, i - 1])

    return fluxes


def specify_boundary_condition(U_initial):
    """Specify BC."""
    F_left = calculate_cell_F(U_initial[:, 0])
    F_right = calculate_cell_F(U_initial[:, -1])
    return F_left, F_right


def advance(U, facefluxes, dt, dx, N, timestep='euler'):
    """Time integration. Finite volume method."""
    if timestep == 'euler':
        for i in range(N):
            U[:, i] = U[:, i] + dt / dx * (facefluxes[:, i] - facefluxes[:, i + 1])

    return U


def plot_U(x, U, **kwargs):
    """Plot the solution."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(nrows=3, ncols=1)

    ax[0].plot(x, U[0, :], **kwargs)
    ax[0].set_ylabel('Rho')
    
    ax[1].plot(x, U[1, :]/U[0, :], **kwargs)
    ax[1].set_ylabel('U')
    
    ax[2].plot(x, (GAMMA - 1) * (U[2, :] - 0.5 * U[1, :]**2 / U[0, :]), **kwargs)
    ax[2].set_ylabel('p')
    ax[2].set_xlabel('x')
    plt.show()

# --------------------------------------------------------------------------------
if __name__ == '__main__':

    N = 400

    # Inputs
    xmin = 0.0
    xmax = 1.0
    tmax = 0.001
    dt = 0.001
    flux = 'Rusanov'
    
    
    # Discretize domain
    x = np.linspace(xmin ,xmax, N+1)  # Position of faces.
    dx = np.abs(x[1] - x[0])          # Cellsize
    
    # Initialize U at t=0
    U = set_initial_state(N, sampletest=4)
    U_initial = U * 1.0
    
    # Main loop
    t = 0.0
    while t < tmax:

        F = np.zeros((3, N + 1))
        for i in range(N):
            F[:, i] = calculate_cell_F(U[:, i])
               
        facefluxes = calculate_fluxes(U, F, N, typeflux=flux)        
        facefluxes[:, 0], facefluxes[:, -1] = specify_boundary_condition(U_initial)
        
        U = advance(U, facefluxes, dt, dx, N, timestep='euler')

        t += dt


    # Plot solution at t=0 and t=tmax
    xc = np.zeros(N)
    for i in range(N):
        xc[i] = 0.5 * (x[i] + x[i + 1])

    plot_U(xc, U_initial, color='black', ls='dotted')
    plot_U(xc, U, color='red')

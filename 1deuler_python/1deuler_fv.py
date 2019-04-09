#!/usr/bin/env python3
"""
Solve the 1D Euler equations using the 
finite-volume method.

So far this is first order because the solution
in each volume is only represented by cell averages.
Higher orders may require the use of a limiter.
"""
import numpy as np

# CONSTANTS
GAMMA = 1.4


def calc_p(u):
    """Given a state vector calculate pressure"""
    return (GAMMA - 1) * (u[2] - 0.5 * u[1]**2 / u[0])


def calc_F(u):
    """Given a a state vector U calculate the flux vector F"""
    F = np.zeros(3)
    p = calc_p(u)
    F[0] = u[1]
    F[1] = u[1]**2 / u[0] + p
    F[2] = u[1] / u[0] * (u[2] + p)
    return F
    
def calc_Fhat(ul, ur, typeflux='hlle', dx=None, dt=None):
    """Calculate the flux at the interface."""
    Fl = calc_F(ul)
    Fr = calc_F(ur)

    if 'rusanov' in typeflux:
        pl = calc_p(ul)
        pr = calc_p(ur)

        slmax = np.maximum(0., -ul[1] / ul[0] + np.sqrt(GAMMA * pl / ul[0]))
        srmax = np.maximum(0., ur[1] / ur[0] + np.sqrt(GAMMA * pr / ur[0]))
        
        smax = np.maximum(slmax, srmax)

        return 0.5 * (Fl + Fr) - 0.5 * smax * (ur - ul)
        
    elif 'roe' in typeflux:

        # Calculate Roe averages
        sqrtrhol = np.sqrt(ul[0])
        sqrtrhor = np.sqrt(ur[0])
        Hr = ur[2] / ur[0] + calc_p(ur) / ur[0] 
        Hl = ul[2] / ul[0] + calc_p(ul) / ul[0]

        ubar = (sqrtrhol * ul[1] / ul[0] + sqrtrhor * ur[1] / ur[0]) / \
               (sqrtrhol + sqrtrhor)

        Hbar = (sqrtrhol * Hl + sqrtrhor * Hr) / \
               (sqrtrhol + sqrtrhor)

        cbar = np.sqrt(GAMMA * (GAMMA - 1) / (2 - GAMMA) * (Hbar - 0.5 * ubar**2))

        R = np.zeros((3, 3))
        R[0, :] = 1
        R[1, 0] = ubar - cbar
        R[1, 1] = ubar
        R[1, 2] = ubar + cbar
        R[2, 0] = Hbar - ubar * cbar
        R[2, 1] = 0.5 * ubar**2
        R[2, 2] = Hbar + ubar * cbar
        L = np.linalg.inv(R)
        Lambda = np.abs(np.diag(np.array([ubar - cbar, ubar, ubar + cbar])))

        # Entropy fix
        epsilon = 0.05 * cbar
        for i in range(3):
            Lambda[i, i] = (epsilon**2 + Lambda[i, i]**2) / (2 * epsilon) \
                           if np.abs(Lambda[i, i]) < epsilon else Lambda[i, i]
            
        return 0.5 * (Fl + Fr) - 0.5 *\
            np.matmul(R, np.matmul(Lambda, np.matmul(L, ur - ul)))           

    elif 'laxfriedrichs' in typeflux:
        return 0.5 * (Fl + Fr) - 0.5 * dx / dt * (ur - ul)

    elif 'hlle' in typeflux:

        pl = calc_p(ul)
        pr = calc_p(ur)
        
        slmax = np.maximum(0., -ul[1] / ul[0] + np.sqrt(GAMMA * pl / ul[0]))
        slmin = np.minimum(0., -ul[1] / ul[0] - np.sqrt(GAMMA * pl / ul[0]))
        srmax = np.maximum(0., ur[1] / ur[0] + np.sqrt(GAMMA * pr / ur[0]))
        srmin = np.minimum(0., ur[1] / ur[0] - np.sqrt(GAMMA * pr / ur[0]))        
        smin = np.minimum(slmin, srmin)
        smax = np.maximum(slmax, srmax)
        
        return 0.5 * (Fl + Fr) - 0.5 * ((smax + smin) / (smax - smin)) * (Fr - Fl) \
               + ((smax * smin) / (smax - smin)) * (ur - ul) 
    
    elif 'linde' in typeflux:
        pass


def advance_in_time_explicit(u, Flface, Frface, dt, dx, tsteptype='euler'):
    """ Time marching (explicit)"""
    if 'euler' in tsteptype:
        return u + dt / dx * (Flface - Frface)


def add_axes(ax, x, u, **kwargs):
    """Make one plot and return the axes"""
    ax.plot(x, u, **kwargs)
    return ax

if __name__ == "__main__":

    # Setup
    N = 400
    xmin = 0.0
    xmax = 1.0
    tmax = 0.25
    dt = 0.001
    x = np.linspace(xmin, xmax, N + 1)
    dx = x[1] - x[0]
    u = np.zeros((3, N))
    Fhat = np.zeros((3, N + 1))
    tflux = 'roe'
    

    # Initial condition
    midway = int(N/2) if N % 2 == 0 else int(N/2 + 1)
    u[0, :midway] = 1.0
    u[1, :midway] = 0.0
    u[2, :midway] = 1.0
    u[0, midway:] = 0.125
    u[1, midway:] = 0.0
    u[2, midway:] = 0.1
    u_init = u * 1.0
    
    # Update states until tmax
    t = 0
    while t < tmax:
        # Boundary condition
        Fhat[:, 0] = calc_F(u_init[:, 0])
        Fhat[:, -1] = calc_F(u_init[:, -1]) 
        for i in range(1, N):
            Fhat[:, i] = calc_Fhat(u[:, i - 1], u[:, i], typeflux=tflux, dx=dx, dt=dt)
            u[:, i - 1] = advance_in_time_explicit(u[:, i - 1], Fhat[:, i - 1], Fhat[:, i], dt, dx)

        t += dt

    # Calculate fine grid solution
    
    # Plot solution
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(3, 1, figsize=(11, 8.5))
    xc = np.zeros(N)
    for i in range(N):
        xc[i] = 0.5 * (x[i] + x[i + 1])

    ax[0].scatter(xc, u[0, :], color='k')
    ax[1].scatter(xc, u[1, :] / u[0, :], color='k')
    ax[2].scatter(xc, calc_p(u), color='k')
    
    for i in range(3):
        ax[i].set_xlabel('x')

    plt.show()
    
    

#!/usr/bin/env python3
"""
Solve the 1D Euler equations using the 
finite-volume method.
"""
import numpy as np
import matplotlib.pyplot as plt

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

        slmax = np.max(0, -ul[1] / ul[0] + np.sqrt(GAMMA * pl / ul[0]))
        srmax = np.max(0, ur[1] / ur[0] + np.sqrt(GAMMA * pr / ur[0]))
        
        smax = np.max(slmax, srmax)

        return 0.5 * (Fl + Fr) - 0.5 * smax * (ur - ul)
        
    elif 'roe' in typeflux:

        # Calculate Roe averages
        ubar = (np.sqrt(ul[0]) * ul[1] / ul[0] + np.sqrt(ur[0]) * ur[1] / ur[0]) / \
               (np.sqrt(ul[0]) + np.sqrt(ur[0]))

        Hr = ur[2] / ur[0]**2 + calc_p(ur) / ur[0] 
        Hl = ul[2] / ul[0]**2 + calc_p(ul) / ul[0]
        
        Hbar = (np.sqrt(ul[0]) * Hl + np.sqrt(ur[0]) * Hr) / \
               (np.sqrt(ul[0]) + np.sqrt(ur[0]))

        cbar = np.sqrt((GAMMA - 1) * (Hbar - 0.5 * ubar**2))

        R = np.zeros((3, 3))
        R[0, :] = 1
        R[1, 0] = ubar - cbar
        R[1, 1] = ubar
        R[1, 2] = ubar + cbar
        R[2, 0] = Hbar - ubar * cbar
        R[2, 1] = 0.5 * ubar**2
        R[2, 2] = Hbar + ubar * cbar
        L = np.linalg.inv(R)
        Lambda = np.diag(np.array([ubar - cbar, ubar, ubar + cbar]))

        return 0.5 * (Fl + Fr) - 0.5 *\
            np.matmul(np.matmul(np.matmul(R, Lambda), L), ur - ul)           

    elif 'laxfriedrichs' in typeflux:
        return 0.5 * (Fl + Fr) - 0.5 * dx / dt * (ur - ul)

    elif 'hlle' in typeflux:
        
        slmax = np.max(0, -ul[1] / ul[0] + np.sqrt(GAMMA * pl / ul[0]))
        slmin = np.min(0, -ul[1] / ul[0] - np.sqrt(GAMMA * pl / ul[0]))
        srmax = np.max(0, ur[1] / ur[0] + np.sqrt(GAMMA * pr / ur[0]))
        srmin = np.min(0, ur[1] / ur[0] - np.sqrt(GAMMA * pr / ur[0]))        
        smin = np.min(slmin, srmin)
        smax = np.max(slmax, srmax)
        
         return 0.5 * (Fl + Fr) - 0.5 * ((smax + smin) / (smax - smin)) * (Fr - Fl) \
               + ((smax * smin) / (smax - smin)) * (ur - ul) 

    elif 'linde' in typeflux:
        pass


def advance_in_time_explicit(u, Flface, Frface, dt, dx, tsteptype='euler'):
    """ Time marching (explicit)"""
    if 'euler' in tsteptype:
        return u + dt / dx * (Flface - Frface)



if __name__ == "__main__":



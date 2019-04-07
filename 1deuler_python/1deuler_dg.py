#!/usr/bin/env python3
"""
Solve the 1D Euler equations using the 
Discontinuous-Galerkin formulation.

Algorithm - 
    1.  
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
        diffusionterm = np.zeros(3)

        
    

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

    elif 'linde' in typeflx:
        pass



    



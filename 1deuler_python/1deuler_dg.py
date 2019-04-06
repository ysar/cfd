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


def calc_F(U):
    """Given a a state vector U calculate the flux vector F"""
    pass


def calc_R():
    """Calculate the residual for the DG formulation."""
    pass


def calc_Fhat(leftstate, rightstate):
    """Calculate the flux at the interface."""
    pass


def calc_phifaces():
    """DG method needs calculation of the basis function at interfaces"""
    pass


def calc_M():
    """Calculate the mass matrix. Size (N*(p+1), s)"""
    pass


def unroll_to_vector(matrix):
    """Given the MU matrix of size (N*(p+1), s), unroll it into a 
    vector of size N*(p+1)*s."""
    pass


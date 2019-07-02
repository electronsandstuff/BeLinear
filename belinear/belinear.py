#!/usr/bin/env python

################################################################################
# File: belinear.py
# Description: belinear provides methods that make working with linear beam
#              optics easy and fast
# Author: Christopher M. Pierce (cmp285@cornell.edu)
################################################################################

################################################################################
# Imports
################################################################################
import numpy as np
import functools as f

################################################################################
# Functions
################################################################################
def get_w_from_gamma(g):
    return np.arcsinh(np.sqrt(g**2 - 1))

def get_delta_w_from_gamma(g, dg):
    return np.log((g+dg+np.sqrt((g+dg)**2-1))/(g+np.sqrt(g**2-1)))

def get_theta_prime(E, B, gamma_initial, delta_z):
    # Set some constants
    c = 299792458 # Speed of light
    mc2 = 511e3 # Rest mass of electrons

    # Convert to normalized fields
    E_norm = E/mc2
    b = B/(2*mc2)*c

    # Get Delta gamma
    delta_gamma = E_norm*delta_z

    # Find the rapidity
    w  = get_w_from_gamma(gamma_initial)
    dw = get_delta_w_from_gamma(gamma_initial, delta_gamma)

    # Find the larmor rotation
    delta_theta  = 1/(np.real(np.sinc(dw/2/np.pi*1j))*np.sinh(dw/2+w))

    # Return it
    return delta_theta

def get_M_const(E, B, gamma_initial, delta_z):
    # Set some constants
    c = 299792458 # Speed of light
    mc2 = 511e3 # Rest mass of electrons

    # Convert to normalized fields
    E_norm = E/mc2
    b = B/(2*mc2)*c

    # Get theta_prime
    theta_prime = get_theta_prime(E, B, gamma_initial, delta_z)

    # Return the matrix
    return np.array([[
        np.cos(b*delta_z*theta_prime),
        np.sinc(b*delta_z*theta_prime/np.pi)*theta_prime*delta_z
    ],[
        -1*b*np.sin(b*delta_z*theta_prime),
        np.cos(b*delta_z*theta_prime)
    ]])

def get_M_edge(E, gamma, side):
    # Set some constants
    c = 299792458 # Speed of light
    mc2 = 511e3 # Rest mass of electrons

    # Convert to normalized fields
    E_norm = E/mc2

    # Get beta
    beta = np.sqrt(1 - 1/gamma**2)

    # Find the sign constant
    c = 0.0
    if(side == 'rising'):
        c = -1.0
    elif(side == 'falling'):
        c = 1.0
    else:
        raise ValueError('Unknown side')

    # Find the matrix
    n = E.shape[0]
    M = np.reshape(np.repeat(np.identity(2), n),  (2,2,n))
    M[1, 0, :] = c*E_norm/2/beta
    return M

def get_gammas(E, B, delta_z, gamma_initial = 1):
    # Set some constants
    c = 299792458 # Speed of light
    mc2 = 511e3 # Rest mass of electrons

    # Convert to normalized fields
    E_norm = E/mc2

    # Add the first section (no change in gamma)
    E_norm = np.insert(E_norm, 0, 0.0)

    # Calculate the gammas
    gammas = np.cumsum(E_norm)*delta_z + gamma_initial

    return gammas

def get_M(E, B, delta_z):
    # Set the initial conditions for the transfer matrix calculation
    gamma = 1.0

    # Get the gammas
    gammas = get_gammas(E, B, delta_z)

    # Calculate the constant matrices and edge matrices
    Ms = get_M_const(E, B, gammas[:-1], delta_z)
    rising_Ms = np.concatenate((np.array([[[1.0,], [0.0,]], [[0.0,], [1.0,]]]), get_M_edge(E[1:], gammas[1:-1], 'rising')), axis=2)
    falling_Ms = get_M_edge(E, gammas[1:], 'falling')

    # Interleave the arrays
    c = np.empty((2,2, Ms.shape[-1]+rising_Ms.shape[-1]+falling_Ms.shape[-1],), dtype=Ms.dtype)
    c[:,:,0::3] = rising_Ms
    c[:,:,1::3] = Ms
    c[:,:,2::3] = falling_Ms

    # Multiply everything
    return f.reduce(np.dot, c.T).T

def get_B_row(E, B, delta_z):
    # Get m
    M = get_M(E, B, delta_z)

    # Return the array
    return np.array([M[0,0]**2, M[0,1]**2])

def get_B(voltages, anode_map_filename, beamline_length, N=5000):
    # Load the field-maps
    anode_pos = np.genfromtxt(anode_map_filename)[1:][:,0]
    anode_fld = np.genfromtxt(anode_map_filename)[1:][:,1]/9000.0

    # Create the positions we will evaluate the transfer matrices at
    z = np.linspace(0.0, beamline_length, N)
    z = (z[:-1] + z[1:])/2 # Change to the centers of the boxes
    delta_z = z[1]-z[0]

    # Get the E field and B field at each location
    E = np.interp(z, anode_pos, anode_fld, left=0.0, right=0.0)
    B = np.array([0.0 for _ in E])

    return np.array([get_B_row(E*v, B, delta_z) for v in voltages])

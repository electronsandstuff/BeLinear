#!/usr/bin/env python

'''
BeLinear - Performant numerical solutions of the paraxial ray equation
Copyright (C) 2020 Christopher M. Pierce (contact@chris-pierce.com)

This program is free software; you can redistribute it and/or
modify it under the terms of version 3 of the GNU Lesser General Public
License as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program; if not, write to the Free Software Foundation,
Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
'''

import numpy as np
import matprod
import scipy.integrate as integrate
from . import tests

################################################################################
# Physical Constants
################################################################################
m = 9.1093837015e-31  # Mass of the electron in kg
c = 299792458         # Speed of light in m/s
q = 1.602176634e-19   # Fundamental charge in C
mc2 = 510.998950e3    # Mass of the electron in eV

################################################################################
# Low level IVP solving routines
################################################################################
def get_dM_implicit_euler(Ez, Bz, delta_z, gamma_initial=1):
    # Calculate longitudinal energy/velocity along the beamline
    # Use fixed spacing for performance reasons
    gamma_prime = Ez/mc2
    gamma_prime_prime = np.gradient(gamma_prime, delta_z)[1:]
    gamma = gamma_initial + integrate.cumtrapz(gamma_prime, None, dx=delta_z,
                                                   initial=0)[1:]
    beta = np.sqrt(1-1/gamma**2)

    # Compute the larmor frequency
    omega_L = (q/2/gamma/m)*Bz[1:]

    # Compute the matrix dM = (1-h*M(z))^(-1) where M(z) is the matrix of
    # coefficients to the coupled ODE in the expression y'(z) = M @ y(z) and h
    # is the step size (delta_z).  The following is for the set of coefficients
    # in the paraxial ray equation with y = [x, P_x] and
    # M = [[0, 1/(beta(z)*gamma(z)) ],
    #      [-(2*gamma(z)*omega_L^2 + c^2*gamma''(z)) / (2c^2*beta(z)) ]]
    denom = 2*gamma*(c**2+delta_z**2/beta**2*omega_L**2)
    denom = denom + c**2*delta_z**2/beta**2*gamma_prime_prime
    dM = np.empty((2,2,Ez.size))
    dM[0,0,1:] = 2*c**2*gamma/denom
    dM[0,1,1:] = 2*c**2*delta_z/beta/denom
    dM[1,0,1:] = -1*delta_z/beta*gamma
    dM[1,0,1:] = dM[1,0,1:]*(2*gamma*omega_L**2 + c**2*gamma_prime_prime)/denom
    dM[0,0,0]  = 1.0
    dM[1,0,0]  = 0.0
    dM[0,1,0]  = 0.0
    dM[1,1,:]  = dM[0,0,:]

    # Return the matrices
    return dM

def get_dM_midpoint(Ez, Bz, delta_z, gamma_initial=1):
    # Calculate longitudinal energy/velocity along the beamline
    # Use fixed spacing for performance reasons
    gamma_prime = Ez/mc2
    gamma_prime_prime = np.gradient(gamma_prime, delta_z)
    gamma = gamma_initial + integrate.cumtrapz(gamma_prime, None, dx=delta_z,
                                                   initial=0)
    beta = np.sqrt(1-1/gamma**2)

    # Compute the larmor frequency
    omega_L = (q/2/gamma/m)*Bz

    # Transform them all to z+1/2 spacing
    gamma_prime_prime = (gamma_prime_prime[1:] + gamma_prime_prime[:-1])/2
    gamma = (gamma[1:] + gamma[:-1])/2
    beta = (beta[1:] + beta[:-1])/2
    omega_L = (omega_L[1:] + omega_L[:-1])/2

    # Compute the matrix T = (1-h/2*M(z))^(-1)(1+h/2*M(z)) where M(z) is the
    # matrix of coefficients to the coupled ODE in the expression
    # y'(z) = M @ y(z) and h is the step size (delta_z).  The following is for
    # the set of coefficients in the paraxial ray equation with y = [x, P_x] and
    # M = [[0, 1/(beta(z)*gamma(z)) ],
    # [-(2*gamma(z)*omega_L^2 + c^2*gamma''(z)) / (2c^2*beta(z)) ]]
    denom = 2*gamma*(4*c**2*beta**2 + delta_z**2*omega_L**2)
    denom = denom + c**2*delta_z**2*gamma_prime_prime
    dM = np.empty((2,2,Ez.size))
    dM[0,0,1:] = 16*c**2*beta**2*gamma/denom - 1
    dM[0,1,1:] = 8*c**2*delta_z*beta/denom
    dM[1,0,1:] = -4*delta_z*beta*gamma
    dM[1,0,1:] = dM[1,0,1:]*(2*gamma*omega_L**2 + c**2*gamma_prime_prime)/denom
    dM[0,0,0]  = 1.0
    dM[1,0,0]  = 0.0
    dM[0,1,0]  = 0.0
    dM[1,1,:]  = dM[0,0,:]

    # Return the matrices
    return dM

def get_dM_constant_field(Ez, Bz, delta_z, gamma_initial=1):
    # Compute the larmor frequency
    b = Bz[:-1]/(2*mc2)*c
    gamma_prime = Ez/mc2

    # Calculate longitudinal energy/velocity along the beamline
    # Use fixed spacing for performance reasons
    gamma = gamma_initial + integrate.cumtrapz(gamma_prime, None, dx=delta_z,
                                                   initial=0)
    gamma_prime_prime = np.gradient(gamma_prime, delta_z)[:-1]
    beta = np.sqrt(1-1/gamma**2)
    w = np.arcsinh(np.sqrt(gamma**2 - 1))

    # Compute the scaling factor for the larmor rotation in each step
    ls = 1/(np.real(np.sinc((w[1:] - w[:-1])/2/np.pi*1j)))
    ls = ls/np.sinh((w[1:] + w[:-1])/2)

    # Calculate the transfer matrices for the constant field regions plus on
    # falling and rising edge.  IE T = M_rise @ M_fall @ M_const.
    COS = np.cos(b*delta_z*ls)
    SIN = -b*np.sin(b*delta_z*ls)
    SNC = np.sinc(b*delta_z*ls/np.pi)*ls*delta_z
    dM = np.empty((2,2,Ez.size))
    dM[0,0,1:] = COS
    dM[0,1,1:] = SNC
    dM[1,0,1:] = SIN - COS*gamma_prime_prime*delta_z/2/beta[1:]
    dM[1,1,1:] = COS - SNC*gamma_prime_prime*delta_z/2/beta[1:]
    dM[0,0,0]  = 1.0
    dM[1,0,0]  = 0.0
    dM[0,1,0]  = 0.0
    dM[1,1,0]  = 1.0

    # Multiply the matrices together
    return dM

def get_dM(Ez, Bz, delta_z, gamma_initial=1, method='midpoint'):
    # Get the matrices dM
    if(method == 'midpoint'):
        dM = get_dM_midpoint(Ez, Bz, delta_z, gamma_initial)
    elif(method == 'implicit_euler'):
        dM = get_dM_implicit_euler(Ez, Bz, delta_z, gamma_initial)
    elif(method == 'constant_field'):
        dM = get_dM_constant_field(Ez, Bz, delta_z, gamma_initial)
    else:
        raise ValueError('Unknown solver method "{:s}"'.format(method))

    # Return it
    return dM

################################################################################
# High level interface
################################################################################
def get_M(z, Ez, Bz, h, gamma_initial=1, method='midpoint'):
    # Interpolate Ez and Bz onto the requested step size
    dz  = np.arange(min(z), max(z), h)
    dEz = np.interp(dz, z, Ez)
    dBz = np.interp(dz, z, Bz)

    # Return the product of the matrices dM
    return matprod.lprod(get_dM(dEz, dBz, h, gamma_initial, method))

def get_cum_M(z, Ez, Bz, h, gamma_initial=1, method='midpoint'):
    # Interpolate Ez and Bz onto the requested step size
    dz  = np.arange(min(z), max(z), h)
    dEz = np.interp(dz, z, Ez)
    dBz = np.interp(dz, z, Bz)

    # Return the cumulative product of the matrices dM
    return matprod.cumlprod(get_dM(dEz, dBz, h, gamma_initial, method))

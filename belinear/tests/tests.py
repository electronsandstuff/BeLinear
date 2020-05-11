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

################################################################################
# Imports
################################################################################
import belinear
import pickle
import numpy as np
import unittest
import os

################################################################################
# The tests
################################################################################
class TestBeLinear(unittest.TestCase):
    def setUp(self):
        '''Loads the data for tests'''
        # Get our path and the path of the data file
        this_path = os.path.dirname(__file__)
        filename = os.path.join(this_path, 'test_data.pickle')

        # Load the test data
        with open(filename, 'rb') as f:
            self.test_data = pickle.load(f)

    def test_get_M_midpoint_rest(self):
        # Pull out the test data
        z   = self.test_data['z_raw']
        Ez  = self.test_data['Ez_raw']
        Bz  = self.test_data['Bz_raw']
        valid_M = self.test_data['M1']

        # Get the matrix
        h = 5e-6
        test_M = belinear.get_M(z, Ez, Bz, h, gamma_initial=1,
                                method='midpoint')

        # Compare the matrices
        self.assertTrue(np.isclose(valid_M, test_M).all())

    def test_get_M_euler_rest(self):
        # Pull out the test data
        z   = self.test_data['z_raw']
        Ez  = self.test_data['Ez_raw']
        Bz  = self.test_data['Bz_raw']
        valid_M = self.test_data['M1']

        # Get the matrix
        h = 5e-6
        test_M = belinear.get_M(z, Ez, Bz, h, gamma_initial=1,
                                method='implicit_euler')

        # Compare the matrices
        self.assertTrue(np.isclose(valid_M, test_M, rtol=0.1).all())

    def test_get_M_constant_rest(self):
        # Pull out the test data
        z   = self.test_data['z_raw']
        Ez  = self.test_data['Ez_raw']
        Bz  = self.test_data['Bz_raw']
        valid_M = self.test_data['M1']

        # Get the matrix
        h = 5e-6
        test_M = belinear.get_M(z, Ez, Bz, h, gamma_initial=1,
                                method='constant_field')

        # Compare the matrices
        self.assertTrue(np.isclose(valid_M, test_M, rtol=0.01).all())

    def test_get_M_midpoint_rel(self):
        # Pull out the test data
        z   = self.test_data['z_raw']
        Ez  = self.test_data['Ez_raw']
        Bz  = self.test_data['Bz_raw']
        valid_M = self.test_data['M2']

        # Get the matrix
        h = 5e-6
        test_M = belinear.get_M(z, Ez, Bz, h, gamma_initial=2,
                                method='midpoint')

        # Compare the matrices
        self.assertTrue(np.isclose(valid_M, test_M).all())

    def test_get_M_euler_rel(self):
        # Pull out the test data
        z   = self.test_data['z_raw']
        Ez  = self.test_data['Ez_raw']
        Bz  = self.test_data['Bz_raw']
        valid_M = self.test_data['M2']

        # Get the matrix
        h = 5e-6
        test_M = belinear.get_M(z, Ez, Bz, h, gamma_initial=2,
                method='implicit_euler')

        # Compare the matrices
        self.assertTrue(np.isclose(valid_M, test_M, rtol=0.1).all())

    def test_get_M_constant_rel(self):
        # Pull out the test data
        z   = self.test_data['z_raw']
        Ez  = self.test_data['Ez_raw']
        Bz  = self.test_data['Bz_raw']
        valid_M = self.test_data['M2']

        # Get the matrix
        h = 5e-6
        test_M = belinear.get_M(z, Ez, Bz, h, gamma_initial=2,
                method='constant_field')

        # Compare the matrices
        self.assertTrue(np.isclose(valid_M, test_M, rtol=0.01).all())

    def test_get_M_cum_rest(self):
        # Pull out the test data
        z   = self.test_data['z_raw']
        Ez  = self.test_data['Ez_raw']
        Bz  = self.test_data['Bz_raw']
        valid_M = self.test_data['M1_cum']

        # Get the matrix
        h = 5e-6
        test_M = belinear.get_cum_M(z, Ez, Bz, h, gamma_initial=1,
                method='midpoint')

        # Compare the matrices
        self.assertTrue(np.isclose(valid_M, test_M).all())

    def test_get_M_cum_rel(self):
        # Pull out the test data
        z   = self.test_data['z_raw']
        Ez  = self.test_data['Ez_raw']
        Bz  = self.test_data['Bz_raw']
        valid_M = self.test_data['M2_cum']

        # Get the matrix
        h = 5e-6
        test_M = belinear.get_cum_M(z, Ez, Bz, h, gamma_initial=2,
                method='midpoint')

        # Compare the matrices
        self.assertTrue(np.isclose(valid_M, test_M).all())

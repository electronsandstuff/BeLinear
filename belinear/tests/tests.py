################################################################################
# Unit tests for belinear- Christopher M. Pierce (chris@chris-pierce.com)
################################################################################

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
        z_raw   = self.test_data['z_raw']
        Ez_raw  = self.test_data['Ez_raw']
        Bz_raw  = self.test_data['Bz_raw']
        valid_M = self.test_data['M1']

        # Interpolate to the points used in ODE solving
        z       = np.linspace(z_raw.min(), z_raw.max(), 100000)
        delta_z = z[1] - z[0]
        Ez      = np.interp(z, z_raw, Ez_raw)
        Bz      = np.interp(z, z_raw, Bz_raw)

        # Get the matrix
        test_M = belinear.get_M(Ez, Bz, delta_z, gamma_initial=1,
                method='midpoint')

        # Compare the matrices
        self.assertTrue(np.isclose(valid_M, test_M).all())

    def test_get_M_euler_rest(self):
        # Pull out the test data
        z_raw   = self.test_data['z_raw']
        Ez_raw  = self.test_data['Ez_raw']
        Bz_raw  = self.test_data['Bz_raw']
        valid_M = self.test_data['M1']

        # Interpolate to the points used in ODE solving
        z       = np.linspace(z_raw.min(), z_raw.max(), 100000)
        delta_z = z[1] - z[0]
        Ez      = np.interp(z, z_raw, Ez_raw)
        Bz      = np.interp(z, z_raw, Bz_raw)

        # Get the matrix
        test_M = belinear.get_M(Ez, Bz, delta_z, gamma_initial=1,
                method='implicit_euler')

        # Compare the matrices
        self.assertTrue(np.isclose(valid_M, test_M, rtol=0.1).all())

    def test_get_M_constant_rest(self):
        # Pull out the test data
        z_raw   = self.test_data['z_raw']
        Ez_raw  = self.test_data['Ez_raw']
        Bz_raw  = self.test_data['Bz_raw']
        valid_M = self.test_data['M1']

        # Interpolate to the points used in ODE solving
        z       = np.linspace(z_raw.min(), z_raw.max(), 100000)
        delta_z = z[1] - z[0]
        Ez      = np.interp(z, z_raw, Ez_raw)
        Bz      = np.interp(z, z_raw, Bz_raw)

        # Get the matrix
        test_M = belinear.get_M(Ez, Bz, delta_z, gamma_initial=1,
                method='constant_field')

        # Compare the matrices
        self.assertTrue(np.isclose(valid_M, test_M, rtol=0.01).all())

    def test_get_M_midpoint_rel(self):
        # Pull out the test data
        z_raw   = self.test_data['z_raw']
        Ez_raw  = self.test_data['Ez_raw']
        Bz_raw  = self.test_data['Bz_raw']
        valid_M = self.test_data['M2']

        # Interpolate to the points used in ODE solving
        z       = np.linspace(z_raw.min(), z_raw.max(), 100000)
        delta_z = z[1] - z[0]
        Ez      = np.interp(z, z_raw, Ez_raw)
        Bz      = np.interp(z, z_raw, Bz_raw)

        # Get the matrix
        test_M = belinear.get_M(Ez, Bz, delta_z, gamma_initial=2,
                method='midpoint')

        # Compare the matrices
        self.assertTrue(np.isclose(valid_M, test_M).all())

    def test_get_M_euler_rel(self):
        # Pull out the test data
        z_raw   = self.test_data['z_raw']
        Ez_raw  = self.test_data['Ez_raw']
        Bz_raw  = self.test_data['Bz_raw']
        valid_M = self.test_data['M2']

        # Interpolate to the points used in ODE solving
        z       = np.linspace(z_raw.min(), z_raw.max(), 100000)
        delta_z = z[1] - z[0]
        Ez      = np.interp(z, z_raw, Ez_raw)
        Bz      = np.interp(z, z_raw, Bz_raw)

        # Get the matrix
        test_M = belinear.get_M(Ez, Bz, delta_z, gamma_initial=2,
                method='implicit_euler')

        # Compare the matrices
        self.assertTrue(np.isclose(valid_M, test_M, rtol=0.1).all())

    def test_get_M_constant_rel(self):
        # Pull out the test data
        z_raw   = self.test_data['z_raw']
        Ez_raw  = self.test_data['Ez_raw']
        Bz_raw  = self.test_data['Bz_raw']
        valid_M = self.test_data['M2']

        # Interpolate to the points used in ODE solving
        z       = np.linspace(z_raw.min(), z_raw.max(), 100000)
        delta_z = z[1] - z[0]
        Ez      = np.interp(z, z_raw, Ez_raw)
        Bz      = np.interp(z, z_raw, Bz_raw)

        # Get the matrix
        test_M = belinear.get_M(Ez, Bz, delta_z, gamma_initial=2,
                method='constant_field')

        # Compare the matrices
        self.assertTrue(np.isclose(valid_M, test_M, rtol=0.01).all())

    def test_get_M_cum_rest(self):
        # Pull out the test data
        z_raw   = self.test_data['z_raw']
        Ez_raw  = self.test_data['Ez_raw']
        Bz_raw  = self.test_data['Bz_raw']
        valid_M = self.test_data['M1_cum']

        # Interpolate to the points used in ODE solving
        z       = np.linspace(z_raw.min(), z_raw.max(), 100000)
        delta_z = z[1] - z[0]
        Ez      = np.interp(z, z_raw, Ez_raw)
        Bz      = np.interp(z, z_raw, Bz_raw)

        # Get the matrix
        test_M = belinear.get_cum_M(Ez, Bz, delta_z, gamma_initial=1,
                method='midpoint')

        # Compare the matrices
        self.assertTrue(np.isclose(valid_M, test_M).all())

    def test_get_M_cum_rel(self):
        # Pull out the test data
        z_raw   = self.test_data['z_raw']
        Ez_raw  = self.test_data['Ez_raw']
        Bz_raw  = self.test_data['Bz_raw']
        valid_M = self.test_data['M2_cum']

        # Interpolate to the points used in ODE solving
        z       = np.linspace(z_raw.min(), z_raw.max(), 100000)
        delta_z = z[1] - z[0]
        Ez      = np.interp(z, z_raw, Ez_raw)
        Bz      = np.interp(z, z_raw, Bz_raw)

        # Get the matrix
        test_M = belinear.get_cum_M(Ez, Bz, delta_z, gamma_initial=2,
                method='midpoint')

        # Compare the matrices
        self.assertTrue(np.isclose(valid_M, test_M).all())
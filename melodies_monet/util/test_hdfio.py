import os
import sys
import logging
import unittest
import numpy as np
import hdfio

class test_hdf(unittest.TestCase):

    def setUp(self):
        self.filename = 'test.hdf'
        self.nx = 2
        self.ny = 3
        self.nz = 5
        self.x = np.linspace(0, 1, self.nx, dtype=np.float64)
        self.y = np.linspace(0, 1, self.ny, dtype=np.float64)
        self.z = np.linspace(0, 1, self.nz, dtype=np.float64)
        xfield, yfield, zfield \
            = np.meshgrid(self.x, self.y, self.z, indexing='ij')
        self.u = np.empty((self.nx, self.ny, self.nz), dtype=np.float32)
        self.v = np.empty((self.nx, self.ny, self.nz), dtype=np.float32)
        self.w = np.empty((self.nx, self.ny, self.nz), dtype=np.float32)
        self.u[:,:,:] = 1
        self.v[:,:,:] = xfield[:,:,:] * yfield[:,:,:] * zfield[:,:,:]
        self.w[:,:,:] = xfield[:,:,:]**2 + yfield[:,:,:]**2

    def test_hdf_create(self):
        fileid = hdfio.hdf_create(self.filename)
        hdfio.hdf_write_coord(fileid, 'x', self.x)
        hdfio.hdf_write_coord(fileid, 'y', self.y)
        hdfio.hdf_write_coord(fileid, 'z', self.z)
        hdfio.hdf_write_field(fileid, 'u', ('x', 'y', 'z'), self.u)
        hdfio.hdf_write_field(fileid, 'v', ('x', 'y', 'z'), self.v)
        hdfio.hdf_write_field(fileid, 'w', ('x', 'y', 'z'), self.w)
        hdfio.hdf_close(fileid)

    def test_hdf_list(self):
        fileid = hdfio.hdf_open(self.filename)
        hdfio.hdf_list(fileid)
        hdfio.hdf_close(fileid)

    def test_hdf_read(self):
        fileid = hdfio.hdf_open(self.filename)
        u = hdfio.hdf_read(fileid, 'u')
        v = hdfio.hdf_read(fileid, 'v')
        w = hdfio.hdf_read(fileid, 'w')
        np.testing.assert_equal(self.u, u)
        np.testing.assert_equal(self.v, v)
        np.testing.assert_equal(self.w, w)
        hdfio.hdf_close(fileid)

    def tearDown(self):
        pass

if __name__ == '__main__':

    logging.basicConfig(
        level=logging.DEBUG, stream=sys.stdout)

    test_suite = unittest.TestSuite([
        unittest.TestLoader().loadTestsFromTestCase(test_hdf)])
    unittest.TextTestRunner(verbosity=2).run(test_suite)

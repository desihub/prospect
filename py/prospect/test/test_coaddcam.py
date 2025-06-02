# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test prospect.coaddcam.
"""
import unittest
from ..coaddcam import index_dichotomy, interp_grid, coadd_brz_cameras, coaddcam_prospect


class TestCoaddcam(unittest.TestCase):
    """Test prospect.coaddcam.
    """

    @classmethod
    def setUpClass(cls):
        pass

    @classmethod
    def tearDownClass(cls):
        pass

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_imports(self):
        """Force pytest to recognize this file as a test module, and
        therefore import the objects above.
        """
        self.assertTrue(callable(index_dichotomy))
        self.assertTrue(callable(interp_grid))
        self.assertTrue(callable(coadd_brz_cameras))
        self.assertTrue(callable(coaddcam_prospect))

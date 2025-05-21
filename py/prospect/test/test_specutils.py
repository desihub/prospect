# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test prospect.specutils.
"""
import unittest
from specutils import SpectrumList
from ..specutils import Spectra, write_spectra, read_spectra, read_spPlate, read_spZbest, read_frame_as_spectra


class TestSpecutils(unittest.TestCase):
    """Test prospect.specutils.
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
        self.assertTrue(issubclass(Spectra, SpectrumList))
        self.assertTrue(callable(write_spectra))
        self.assertTrue(callable(read_spectra))
        self.assertTrue(callable(read_spPlate))
        self.assertTrue(callable(write_spectra))
        self.assertTrue(callable(read_spZbest))
        self.assertTrue(callable(read_frame_as_spectra))

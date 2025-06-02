# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test prospect.viewer.
"""
import unittest
from ..viewer import create_model, plotspectra


class TestViewer(unittest.TestCase):
    """Test prospect.viewer.
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
        self.assertTrue(callable(create_model))
        self.assertTrue(callable(plotspectra))

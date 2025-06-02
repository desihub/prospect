# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test prospect.grid_thumbs.
"""
import unittest
from ..grid_thumbs import grid_thumbs


class TestGridThumbs(unittest.TestCase):
    """Test prospect.grid_thumbs.
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
        self.assertTrue(callable(grid_thumbs))

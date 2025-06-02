# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test prospect.viewer.plots.
"""
import unittest
from ..viewer.plots import _cross_hair_points, _viewer_urls, ViewerPlots


class TestViewerPlots(unittest.TestCase):
    """Test prospect.viewer.plots.
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
        self.assertTrue(issubclass(ViewerPlots, object))
        self.assertTrue(callable(_cross_hair_points))
        self.assertTrue(callable(_viewer_urls))

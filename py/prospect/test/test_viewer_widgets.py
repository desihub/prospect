# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test prospect.viewer.widgets.
"""
import unittest
from ..viewer.widgets import _metadata_table, ViewerWidgets


class TestViewerWidgets(unittest.TestCase):
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
        self.assertTrue(issubclass(ViewerWidgets, object))
        self.assertTrue(callable(_metadata_table))

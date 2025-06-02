# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test prospect.viewer.vi_widgets.
"""
import unittest
from ..viewer.vi_widgets import ViewerVIWidgets


class TestViewerVIWidgets(unittest.TestCase):
    """Test prospect.viewer.vi_widgets.
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
        self.assertTrue(issubclass(ViewerVIWidgets, object))

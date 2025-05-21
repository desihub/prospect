# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test prospect.viewer.cds.
"""
import unittest
from ..viewer.cds import _airtovac, ViewerCDS


class TestViewerCDS(unittest.TestCase):
    """Test prospect.viewer.cds.
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
        self.assertTrue(issubclass(ViewerCDS, object))
        self.assertTrue(callable(_airtovac))

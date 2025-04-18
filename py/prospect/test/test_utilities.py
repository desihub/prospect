# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test prospect.utilities.
"""
import unittest
import re
import sys
from ..utilities import vi_file_fields, get_resources


class TestUtilities(unittest.TestCase):
    """Test prospect.utilities.
    """

    @classmethod
    def setUpClass(cls):
        cls.vi_colnames = [v[0] for v in vi_file_fields]
        cls.vi_dtype = [v[2] for v in vi_file_fields]

    @classmethod
    def tearDownClass(cls):
        pass

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_get_resources(self):
        """Test caching of resource files.
        """
        foo = get_resources('templates')
        self.assertIn('template_index.html', foo.keys())
        self.assertIsInstance(foo['template_index.html'], str)

        bar = get_resources('js')
        self.assertIn('FileSaver.js', bar.keys())
        self.assertIsInstance(bar['FileSaver.js'], str)

        with self.assertRaises(ValueError):
            bad = get_resources('foo')



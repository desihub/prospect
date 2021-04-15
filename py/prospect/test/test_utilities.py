# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test prospect.utilities.
"""
import unittest
import re
import sys
from pkg_resources import resource_filename
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
        bar = get_resources('js')
        with self.assertRaises(ValueError):
            bad = get_resources('foo')


def test_suite():
    """Allows testing of only this module with the command::

        python setup.py test -m <modulename>
    """
    return unittest.defaultTestLoader.loadTestsFromName(__name__)

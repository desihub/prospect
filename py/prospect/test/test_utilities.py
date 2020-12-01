# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test prospect.utilities.
"""
import unittest
import re
import sys
from ..utilities import get_resources


class TestUtilities(unittest.TestCase):
    """Test prospect.utilities.
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

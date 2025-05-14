# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test top-level prospect functions.
"""
import unittest
import re
from .. import __version__ as theVersion


class TestTopLevel(unittest.TestCase):
    """Test top-level prospect functions.
    """

    @classmethod
    def setUpClass(cls):
        cls.versionre = re.compile(
                r'([0-9]+!)?([0-9]+)(\.[0-9]+)*((a|b|rc|\.post|\.dev)[0-9]+)?')

    @classmethod
    def tearDownClass(cls):
        pass

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_version(self):
        """Ensure the version conforms to PEP386/PEP440.
        """
        self.assertRegex(theVersion, self.versionre)



# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test prospect.scripts.
"""
import unittest

_desi_missing = False
try:
    from ..scripts.prospect_pages import _parse, _filter_list, load_spectra_zcat_from_dbentry, page_subset, main
    from ..scripts.prospect_std_templates import main as std_templates_main
except ImportError:
    _desi_missing = True


class TestScripts(unittest.TestCase):
    """Test prospect.scripts.
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

    @unittest.skipIf(_desi_missing, "Skipping test that require DESI software installation.")
    def test_imports(self):
        """Force pytest to recognize this file as a test module, and
        therefore import the objects above.
        """
        self.assertTrue(callable(_parse))
        self.assertTrue(callable(_filter_list))
        self.assertTrue(callable(load_spectra_zcat_from_dbentry))
        self.assertTrue(callable(page_subset))
        self.assertTrue(callable(main))
        self.assertTrue(callable(std_templates_main))

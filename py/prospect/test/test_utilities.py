# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test prospect.utilities.
"""
import unittest
from ..utilities import vi_file_fields, get_resources, get_subset_label


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

    def test_get_subset_label(self):
        """Test subset labels.
        """
        with self.assertRaises(ValueError):
            bad = get_subset_label('20250514', 'unknown')
        self.assertEqual(get_subset_label('20250514', 'cumulative'), 'thru20250514')
        self.assertEqual(get_subset_label('20250514', 'perexp'), 'exp20250514')
        self.assertEqual(get_subset_label('20250514', 'pernight'), '20250514')
        self.assertEqual(get_subset_label('20250514', 'exposures'), '20250514')
        self.assertEqual(get_subset_label('20250514', 'healpix'), '20250514')

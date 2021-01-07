#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
#
# Standard imports
#
import glob
import os
import sys
#
# setuptools' sdist command ignores MANIFEST.in
#
from distutils.command.sdist import sdist as DistutilsSdist
from setuptools import setup, find_packages
#
# DESI support code.
#
have_desiutil = True
try:
    import desiutil.setup as ds
except ImportError:
    have_desiutil = False
#
# Begin setup
#
setup_keywords = dict()
#
# THESE SETTINGS NEED TO BE CHANGED FOR EVERY PRODUCT.
#
setup_keywords['name'] = 'prospect'
setup_keywords['description'] = 'DESI spectrum visualization package'
setup_keywords['author'] = 'DESI Collaboration'
setup_keywords['author_email'] = 'desi-data@desi.lbl.gov'
setup_keywords['license'] = 'BSD'
setup_keywords['url'] = 'https://github.com/desihub/prospect'
#
# END OF SETTINGS THAT NEED TO BE CHANGED.
#
if have_desiutil:
    setup_keywords['version'] = ds.get_version(setup_keywords['name'])
else:
    try:
        with open(os.path.join('py', setup_keywords['name'], '_version.py')) as v:
            setup_keywords['version'] = v.read().split('=')[1].strip().strip("'").strip('"')
    except FileNotFoundError:
        setup_keywords['version'] = '0.0.1'
#
# Use README.rst as long_description.
#
setup_keywords['long_description'] = ''
if os.path.exists('README.rst'):
    with open('README.rst') as readme:
        setup_keywords['long_description'] = readme.read()
#
# Set other keywords for the setup function.  These are automated, & should
# be left alone unless you are an expert.
#
# Treat everything in bin/ except *.rst as a script to be installed.
#
if os.path.isdir('bin'):
    setup_keywords['scripts'] = [fname for fname in glob.glob(os.path.join('bin', '*'))
        if not os.path.basename(fname).endswith('.rst')]
setup_keywords['provides'] = [setup_keywords['name']]
setup_keywords['python_requires'] = '>=3.5'
setup_keywords['zip_safe'] = False
setup_keywords['use_2to3'] = False
setup_keywords['packages'] = find_packages('py')
setup_keywords['package_dir'] = {'': 'py'}
setup_keywords['cmdclass'] = {'sdist': DistutilsSdist}
if have_desiutil:
    setup_keywords['cmdclass']['module_file'] = ds.DesiModule
    setup_keywords['cmdclass']['version'] = ds.DesiVersion
    setup_keywords['cmdclass']['test'] = ds.DesiTest
    setup_keywords['cmdclass']['api'] = ds.DesiAPI
setup_keywords['test_suite']='{name}.test.{name}_test_suite'.format(**setup_keywords)
#
# Autogenerate command-line scripts.
#
# setup_keywords['entry_points'] = {'console_scripts':['run_cmx_htmlfiles = prospect.scripts.prepare_cmx_htmlfiles:main',
#                                                      'run_htmlfiles = prospect.scripts.prepare_htmlfiles:main',
#                                                      'run_specview_cmx_coadds = prospect.scripts.specview_cmx_coadds:main',
#                                                      'run_specview_cmx_frames = prospect.scripts.specview_cmx_frames:main',
#                                                      'run_specview_cmx_targets = prospect.scripts.specview_cmx_targets:main',
#                                                      'run_specview_per_night = prospect.scripts.specview_per_night:main',
#                                                      'run_specview_per_pixel = prospect.scripts.specview_per_pixel:main']}
#
# Add internal data directories.
#
setup_keywords['package_data'] = {'prospect': ['data/*', 'js/*', 'templates/*'],
                                  'prospect.test': ['t/*']}
#
# Run setup command.
#
setup(**setup_keywords)

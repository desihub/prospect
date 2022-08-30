# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
=======================================
prospect.scripts.prospect_std_templates
=======================================

Script to produce "standard templates",
used in prospect.viewer.load_std_templates
"""

import os
import numpy as np
from astropy.table import Table

from desispec.interpolation import resample_flux
from prospect.utilities import load_redrock_templates

#- This script was run once, to produce the following file:
std_template_file = os.path.join(os.environ['HOME'], 'prospect/py/prospect/data/std_templates.fits')
if os.path.isfile(std_template_file):
    print('Error std template file already exists')
    
#- Templates produced from 1st component of old (pre-Aug 2022) Redrock templates:
template_dir = os.path.join(os.environ['DESICONDA'], '../code/redrock-templates/0.7.2')
std_templates = {'QSO': ('QSO',''), 'GALAXY': ('GALAXY',''), 'STAR': ('STAR','F') }
delta_lambd_templates = 3

rr_templts = load_redrock_templates(template_dir=template_dir)
for key,rr_key in std_templates.items() :
    wave_array = np.arange(rr_templts[rr_key].wave[0], rr_templts[rr_key].wave[-1], delta_lambd_templates)
    flux_array = resample_flux(wave_array, rr_templts[rr_key].wave, rr_templts[rr_key].flux[0,:])
    table_templates = Table(data=[wave_array, flux_array], names=['wave_'+key, 'flux_'+key], meta={'name':key})
    table_templates.write(std_template_file, append=True)


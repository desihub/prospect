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
from astropy.io import fits

from desispec.interpolation import resample_flux
from prospect.utilities import load_redrock_templates

# Log of previous versions:
#  std_templates_v0.fits : old Redrock templates only

std_template_file = os.path.join(os.environ['HOME'], 'prospect/py/prospect/data/std_templates.fits')
if os.path.isfile(std_template_file):
    print('Error std template file already exists')
    
#- Templates produced from 1st component of old (pre-Aug 2022) Redrock templates:
template_dir = os.path.join(os.environ['DESICONDA'], '../code/redrock-templates/0.7.2')
#std_templates = {'QSO': ('QSO',''), 'GALAXY': ('GALAXY',''), 'STAR': ('STAR','F') }
std_templates = {'GALAXY': ('GALAXY',''), 'STAR': ('STAR','F') }
delta_lambd_templates = 3

rr_templts = load_redrock_templates(template_dir=template_dir)
for key,rr_key in std_templates.items() :
    wave_array = np.arange(rr_templts[rr_key].wave[0], rr_templts[rr_key].wave[-1], delta_lambd_templates)
    flux_array = resample_flux(wave_array, rr_templts[rr_key].wave, rr_templts[rr_key].flux[0,:])
    table_templates = Table(data=[wave_array, flux_array], names=['wave_'+key, 'flux_'+key], meta={'name':key})
    table_templates.write(std_template_file, append=True)

#- Case of QSO (Summer 2022): use new template provided by A. Brodzeller
qsotemplate_file = os.environ['HOME'] + '/stdtemplate-qso.fits'
hdul = fits.open(qsotemplate_file)
qsowave = 10**(hdul[0].header['CRVAL1']+np.arange(hdul[0].header['NAXIS1'])*hdul[0].header['CDELT1'])
qsoflux = hdul[0].data
resamp_factor = 3
wave_array = qsowave[::resamp_factor]
flux_array = resample_flux(wave_array, qsowave, qsoflux)
table_templates = Table(data=[wave_array, flux_array], names=['wave_QSO', 'flux_QSO'], meta={'name':'QSO'})
table_templates.write(std_template_file, append=True)


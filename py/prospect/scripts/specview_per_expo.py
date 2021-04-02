# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
==================================
prospect.scripts.specview_per_expo
==================================

Script to write html files based on exposures.
It should become useless once specview_per_night works, making use of tiles dir structure.
"""


import os, sys, glob
import argparse
import numpy as np
from astropy.table import Table
import astropy.io.fits

import desispec.io
from desitarget.targetmask import desi_mask
import desispec.spectra
import desispec.frame

from ..myspecselect import myspecselect  # special (to be edited)
from ..viewer import plotspectra
from ..utilities import match_zcat_to_spectra

def _parse():
    parser = argparse.ArgumentParser(description='Create html pages for the spectral viewer')
    parser.add_argument('--datadir', help='Path to spectra relative to DESI_ROOT', type=str, default="datachallenge/reference_runs/18.6/spectro/redux/mini/spectra-64")
    parser.add_argument('--nspecperfile', help='Number of spectra in each html page', type=int, default=50)
    parser.add_argument('--webdir', help='Base directory for webapges', type=str, default=os.environ['HOME']+"/prospect/website")
    parser.add_argument('--vignette_smoothing', help='Smoothing of the vignette images (-1 : no smoothing)', type=float, default=10)
    args = parser.parse_args()

def main():
    args = _parse()
    datadir = os.environ['DESI_ROOT']+"/"+args.datadir
    specfiles = glob.glob(datadir+"/*/*/spectra*")

    specfiles=specfiles[0:2] # TMP

    # Get list of exposures with associated files
    # Todo : select exposures not yet processed
    dict_exposures = {}
    for thespecfile in specfiles :
        spectra = desispec.io.read_spectra(thespecfile)
        pp=np.unique(spectra.fibermap['EXPID'])
        for the_expo in pp :
            if the_expo not in dict_exposures.keys() : dict_exposures[the_expo] = [thespecfile]
            else : dict_exposures[the_expo].append(thespecfile)

    # Loop on exposures
    for exposure,thespecfiles in dict_exposures.items() :
        print("* Working on exposure "+str(exposure))
        print("*   ( Nb of files : "+str(len(thespecfiles))+" )")
        for ii,thefile in enumerate(thespecfiles) :
            thespec = desispec.io.read_spectra(thefile)
            zbfile = thefile.replace('spectra-64-', 'zbest-64-')
            thezb = Table.read(zbfile, 'ZBEST')
            # First creation of spectra with metadata + zbest
            if ii==0 :
                spectra = desispec.spectra.Spectra(meta=thespec.meta, extra=thespec.extra)
                zbest = Table(dtype=thezb.dtype)
            thespec = myspecselect(thespec, expids=[exposure])
            spectra.update(thespec)
            for i_zb in range(len(thezb)) :
                zbest.add_row(thezb[i_zb])

        # Handle several html pages per exposure : sort by fibers
        nbpages = int(np.ceil((spectra.num_spectra()/args.nspecperfile)))
        fiberlist = np.unique(spectra.fibermap["FIBER"])
        if len(fiberlist)!=spectra.num_spectra() : print("!! Several times the same fiber in exposure ??")
        for i_page in range(1,1+nbpages) :
            print("** Page "+str(i_page)+" / "+str(nbpages))
            thespec = myspecselect(spectra, fibers=fiberlist[(i_page-1)*args.nspecperfile:i_page*args.nspecperfile])
            thezb = match_zbest_to_spectra(zbest,thespec)
            titlepage = "specviewer_expo"+str(exposure)+"_fiberset"+str(i_page)
            savedir=args.webdir+"/exposures/expo"+str(exposure)
            if not os.path.exists(savedir) :
                os.mkdir(savedir)
                os.mkdir(savedir+"/vignettes")
            plotspectra(thespec, zcatalog=thezb, title=titlepage, savedir=savedir)
            # for i_spec in range(thespec.num_spectra()) :
            #     saveplot = savedir+"/vignettes/expo"+str(exposure)+"_fiberset"+str(i_page)+"_"+str(i_spec)+".png"
            #     miniplot_spectrum(thespec,i_spec,model=model,saveplot=saveplot, smoothing = args.vignette_smoothing)

# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
===================================
prospect.scripts.specview_per_night
===================================

Write static html files from single-exposure spectra, sorted by night/exposure.
This script uses the prototype tiles datastructure in specprod_dir/tiles.
"""

import os, sys, glob
import argparse
import numpy as np
from astropy.table import Table
import astropy.io.fits

import desispec.io
from desiutil.log import get_logger
from desitarget.targetmask import desi_mask
import desispec.spectra
import desispec.frame

from ..myspecselect import myspecselect # special (to be edited)
from ..viewer import plotspectra
from ..utilities import match_catalog_to_spectra

def _parse():

    parser = argparse.ArgumentParser(description='Create night-based html pages for the spectral viewer')
    parser.add_argument('--specprod_dir', help='overrides $DESI_SPECTRO_REDUX/$SPECPROD/', type=str, default=None)
    parser.add_argument('--nspecperfile', help='Number of spectra in each html page', type=int, default=50)
    parser.add_argument('--webdir', help='Base directory for webapges', type=str, default=None)
    parser.add_argument('--vignette_smoothing', help='Smoothing of the vignette images (-1 : no smoothing)', type=float, default=10)
    args = parser.parse_args()
    return args


def main():
    args = _parse()
    log = get_logger()
    specprod_dir = args.specprod_dir
    if specprod_dir is None : specprod_dir = desispec.io.specprod_root()
    webdir = args.webdir
    if webdir is None : webdir = os.environ["DESI_WWW"]+"/users/armengau/svdc2019c" # TMP, for test

    nights = desispec.io.get_nights(specprod_dir=specprod_dir)
    # TODO - Select night (eg. only last night)
    nights = nights[8:11] # TMP, for test

    for thenight in nights :
        # Get spectra from tiles dir - To consolidate
        specfiles = glob.glob( os.path.join(specprod_dir,"tiles/*/tilespectra-*-"+thenight+".fits") )

        for f in specfiles :
            log.info("Working on file "+f)
            file_label = f[f.find("tilespectra-")+12:f.find(thenight)-1] # From tile-based file description - To consolidate
            spectra = desispec.io.read_spectra(f)
            zbfile = f.replace("tilespectra","zbest")
            zbest = Table.read(zbfile, 'ZBEST')
            # Handle several html pages per pixel : sort by TARGETID
            # NOTE : this way, individual spectra from the same target are together
            # Does it make sense ? (they have the same fit)
            nbpages = int(np.ceil((spectra.num_spectra()/args.nspecperfile)))
            sort_indices = np.argsort(spectra.fibermap["TARGETID"], kind='mergesort') # keep order of equal elts

            for i_page in range(1,1+nbpages) :

                log.info(" * Page "+str(i_page)+" / "+str(nbpages))
                the_indices = sort_indices[(i_page-1)*args.nspecperfile:i_page*args.nspecperfile]
                thespec = myspecselect(spectra, indices=the_indices)
                thezb = match_catalog_to_spectra(zbest,thespec)
                titlepage = "specviewer_night"+thenight+"_"+file_label+"_"+str(i_page)
                html_dir=webdir+"/nights/night"+thenight
                if not os.path.exists(html_dir) :
                    os.makedirs(html_dir)
                    os.mkdir(html_dir+"/vignettes")

                plotspectra(thespec, zcatalog=thezb, title=titlepage, html_dir=html_dir)
#                 for i_spec in range(thespec.num_spectra()) :
#                     saveplot = html_dir+"/vignettes/night"+thenight+"_"+file_label+"_"+str(i_page)+"_"+str(i_spec)+".png"
#                     miniplot_spectrum(thespec, i_spec, model=model, saveplot=saveplot, smoothing = args.vignette_smoothing)

# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
===================================
prospect.scripts.specview_per_pixel
===================================

Write static html files from coadded spectra, sorted by healpixels.
"""

import os, sys, glob
import argparse
import numpy as np
import random
from astropy.table import Table
import astropy.io.fits

import desispec.io
from desiutil.log import get_logger
from desitarget.targetmask import desi_mask
from desitarget.cmx.cmx_targetmask import cmx_mask
from desitarget.sv1.sv1_targetmask import desi_mask as sv1_desi_mask
import desispec.spectra
import desispec.frame

from ..myspecselect import myspecselect  # special (to be edited)
from ..plotframes import coadd_targets, create_model, plotspectra
from ..utilities import specviewer_selection, match_zcat_to_spectra, match_vi_targets  #, miniplot_spectrum


def _parse():

    parser = argparse.ArgumentParser(description='Create pixel-based static html pages for the spectral viewer')
    parser.add_argument('--specprod_dir', help='overrides $DESI_SPECTRO_REDUX/$SPECPROD/', type=str, default=None)
    parser.add_argument('--pixel_list', help='ASCII file providing list of pixels', type=str, default=None)
    parser.add_argument('--mask', help='Select only objects with a given DESI target mask', type=str, default=None)
    parser.add_argument('--gcut', help='Select only objects in a given [dereddened] g-mag range (eg --gcut 22 22.5)', nargs='+', type=float, default=None)
    parser.add_argument('--rcut', help='Select only objects in a given [dereddened] r-mag range (eg --rcut 18 19.5)', nargs='+', type=float, default=None)
    parser.add_argument('--chi2cut', help='Select only objects with Delta_chi2 (from pipeline fit) in a given range (eg --chi2cut 40 100)', nargs='+', type=float, default=None)
    parser.add_argument('--nspecperfile', help='Number of spectra in each html page', type=int, default=50)
    parser.add_argument('--webdir', help='Base directory for webpages', type=str, default=None)
    parser.add_argument('--vignette_smoothing', help='Smoothing of the vignette images (-1 : no smoothing)', type=float, default=10)
    parser.add_argument('--mask_type', help='Mask category : DESI_TARGET,SV1_DESI_TARGET,CMX_TARGET', type=str, default='DESI_TARGET')
    parser.add_argument('--random_pixels', help='Process pixels in random order', action='store_true')
    parser.add_argument('--nmax_spectra', help='Stop the production of HTML pages once a given number of spectra are done', type=int, default=None)
    args = parser.parse_args()
    return args


def main():
    args = _parse()
    log = get_logger()
    specprod_dir = args.specprod_dir
    if specprod_dir is None : specprod_dir = desispec.io.specprod_root()
    webdir = args.webdir
    if webdir is None : webdir = os.environ["DESI_WWW"]+"/users/armengau/svdc2019c" # TMP, for test

    # TODO - Selection on pixels based on existing specviewer pages
    if args.pixel_list is None :
        pixels = glob.glob( os.path.join(specprod_dir,"spectra-64/*/*") )
        pixels = [x[x.rfind("/")+1:] for x in pixels]
    else :
        pixels = np.loadtxt(args.pixel_list, dtype=str)
    if args.random_pixels :
        random.shuffle(pixels)

    # Loop on pixels
    nspec_done = 0
    for pixel in pixels :

        log.info("Working on pixel "+pixel)
        thefile = desispec.io.findfile('spectra', groupname=int(pixel), specprod_dir=specprod_dir)
        individual_spectra = desispec.io.read_spectra(thefile)
        spectra = coadd_targets(individual_spectra)
        zbfile = thefile.replace('spectra-64-', 'zbest-64-')
        if os.path.isfile(zbfile) :
            zbest = Table.read(zbfile, 'ZBEST')
        else :
            log.info("No associated zbest file found : skipping pixel")
            continue

        spectra = specviewer_selection(spectra, log=log,
                        mask=args.mask, mask_type=args.mask_type, gmag_cut=args.gcut, rmag_cut=args.rcut,
                        chi2cut=args.chi2cut, zbest=zbest)
        if spectra == 0 : continue

        # Handle several html pages per pixel : sort by TARGETID
        # TODO - Find a more useful sort ?
        nbpages = int(np.ceil((spectra.num_spectra()/args.nspecperfile)))
        sort_indices = np.argsort(spectra.fibermap["TARGETID"])

        for i_page in range(1,1+nbpages) :

            log.info(" * Page "+str(i_page)+" / "+str(nbpages))
            the_indices = sort_indices[(i_page-1)*args.nspecperfile:i_page*args.nspecperfile]
            thespec = myspecselect(spectra, indices=the_indices)
            thezb, kk = match_zcat_to_spectra(zbest,thespec)
            ### No VI results to display by default
            # VI "catalog" - location to define later ..
            # vifile = os.environ['HOME']+"/prospect/vilist_prototype.fits"
            # vidata = match_vi_targets(vifile, thespec.fibermap["TARGETID"])
            titlepage = "pix"+pixel+"_"+str(i_page)
            if args.gcut is not None :
                titlepage = "gcut-"+str(args.gcut[0])+"-"+str(args.gcut[1])+"_"+titlepage
            if args.rcut is not None :
                titlepage = "rcut-"+str(args.rcut[0])+"-"+str(args.rcut[1])+"_"+titlepage
            if args.chi2cut is not None :
                titlepage = "chi2cut-"+str(args.chi2cut[0])+"-"+str(args.chi2cut[1])+"_"+titlepage
            if args.mask is not None :
                titlepage = args.mask+"_"+titlepage
            model = create_model(thespec, thezb)
            html_dir = os.path.join(webdir,"pix"+pixel)
            if not os.path.exists(html_dir) :
                os.makedirs(html_dir)
                os.mkdir(html_dir+"/vignettes")

            plotspectra(thespec, zcatalog=zbest, model_from_zcat=True, vidata=None, model=None, title=titlepage, html_dir=html_dir, is_coadded=True, mask_type=args.mask_type)
#             for i_spec in range(thespec.num_spectra()) :
#                 saveplot = html_dir+"/vignettes/pix"+pixel+"_"+str(i_page)+"_"+str(i_spec)+".png"
#                 miniplot_spectrum(thespec, i_spec, model=model, saveplot=saveplot, smoothing = args.vignette_smoothing)
            nspec_done += thespec.num_spectra()

        # Stop running if needed, only once a full pixel is completed
        if args.nmax_spectra is not None :
            if nspec_done >= args.nmax_spectra :
                log.info(str(nspec_done)+" spectra done : no other pixel will be processed")
                break

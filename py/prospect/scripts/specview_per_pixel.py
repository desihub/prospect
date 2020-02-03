"""
prospect.scripts.specview_per_pixel
===================================

Write static html files from coadded spectra, sorted by healpixels
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

from prospect import myspecselect # special (to be edited)
from prospect import plotframes
from prospect import utils_specviewer

def parse() :

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


def main(args) :
    
    log = get_logger()
    specprod_dir = args.specprod_dir
    if specprod_dir is None : specprod_dir = desispec.io.specprod_root()
    webdir = args.webdir
    if webdir is None : webdir = os.environ["DESI_WWW"]+"/users/armengau/svdc2019c" # TMP, for test
    if args.mask is not None :
        assert args.mask_type in ['SV1_DESI_TARGET', 'DESI_TARGET', 'CMX_TARGET']
        if args.mask_type ==  'SV1_DESI_TARGET' :
            assert ( args.mask in sv1_desi_mask.names() )
        elif args.mask_type == 'DESI_TARGET' :
            assert ( args.mask in desi_mask.names() ) 
        elif args.mask_type == 'CMX_TARGET' :
            assert ( args.mask in cmx_mask.names() )
            
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
        spectra = plotframes.coadd_targets(individual_spectra)
        zbfile = thefile.replace('spectra-64-', 'zbest-64-')
        if os.path.isfile(zbfile) :
            zbest = Table.read(zbfile, 'ZBEST')
        else :
            log.info("No associated zbest file found : skipping pixel")
            continue

        # Target mask selection
        if args.mask is not None :
            if args.mask_type == 'SV1_DESI_TARGET' :
                w, = np.where( (spectra.fibermap['SV1_DESI_TARGET'] & sv1_desi_mask[args.mask]) )
            elif args.mask_type == 'DESI_TARGET' :
                w, = np.where( (spectra.fibermap['DESI_TARGET'] & desi_mask[args.mask]) )
            elif args.mask_type == 'CMX_TARGET' :
                w, = np.where( (spectra.fibermap['CMX_TARGET'] & cmx_mask[args.mask]) )                
            if len(w) == 0 :
                log.info(" * No "+args.mask+" target in this pixel")
                continue
            else :
                targetids = spectra.fibermap['TARGETID'][w]
                spectra = spectra.select(targets=targetids)

        # Photometry selection
        if args.gcut is not None :
            assert len(args.gcut)==2 # Require range [gmin, gmax]
            gmag = np.zeros(spectra.num_spectra())
            w, = np.where( (spectra.fibermap['FLUX_G']>0) & (spectra.fibermap['MW_TRANSMISSION_G']>0) )
            gmag[w] = -2.5*np.log10(spectra.fibermap['FLUX_G'][w]/spectra.fibermap['MW_TRANSMISSION_G'][w])+22.5
            w, = np.where( (gmag>args.gcut[0]) & (gmag<args.gcut[1]) )
            if len(w) == 0 :
                log.info(" * No target in this pixel with g_mag in requested range")
                continue
            else :
                targetids = spectra.fibermap['TARGETID'][w]
                spectra = spectra.select(targets=targetids)
        if args.rcut is not None :
            assert len(args.rcut)==2 # Require range [rmin, rmax]
            rmag = np.zeros(spectra.num_spectra())
            w, = np.where( (spectra.fibermap['FLUX_R']>0) & (spectra.fibermap['MW_TRANSMISSION_R']>0) )
            rmag[w] = -2.5*np.log10(spectra.fibermap['FLUX_R'][w]/spectra.fibermap['MW_TRANSMISSION_R'][w])+22.5
            w, = np.where( (rmag>args.rcut[0]) & (rmag<args.rcut[1]) )
            if len(w) == 0 :
                log.info(" * No target in this pixel with r_mag in requested range")
                continue
            else :
                targetids = spectra.fibermap['TARGETID'][w]
                spectra = spectra.select(targets=targetids)

        # Chi2 selection
        if args.chi2cut is not None :
            assert len(args.chi2cut)==2 # Require range [chi2min, chi2max]
            thezb, kk = utils_specviewer.match_zcat_to_spectra(zbest,spectra)
            w, = np.where( (thezb['DELTACHI2']>args.chi2cut[0]) & (thezb['DELTACHI2']<args.chi2cut[1]) )
            if len(w) == 0 :
                log.info(" * No target in this pixel with DeltaChi2 in requested range")
                continue
            else :
                targetids = spectra.fibermap['TARGETID'][w]
                spectra = spectra.select(targets=targetids)

        # Handle several html pages per pixel : sort by TARGETID
        # TODO - Find a more useful sort ?
        nbpages = int(np.ceil((spectra.num_spectra()/args.nspecperfile)))
        sort_indices = np.argsort(spectra.fibermap["TARGETID"])
        
        for i_page in range(1,1+nbpages) :
            
            log.info(" * Page "+str(i_page)+" / "+str(nbpages))
            the_indices = sort_indices[(i_page-1)*args.nspecperfile:i_page*args.nspecperfile]
            thespec = myspecselect.myspecselect(spectra, indices=the_indices)
            thezb, kk = utils_specviewer.match_zcat_to_spectra(zbest,thespec)
            ### No VI results to display by default
            # VI "catalog" - location to define later ..
            # vifile = os.environ['HOME']+"/prospect/vilist_prototype.fits"
            # vidata = utils_specviewer.match_vi_targets(vifile, thespec.fibermap["TARGETID"])
            titlepage = "pix"+pixel+"_"+str(i_page)
            if args.gcut is not None :
                titlepage = "gcut-"+str(args.gcut[0])+"-"+str(args.gcut[1])+"_"+titlepage
            if args.rcut is not None :
                titlepage = "rcut-"+str(args.rcut[0])+"-"+str(args.rcut[1])+"_"+titlepage
            if args.chi2cut is not None :
                titlepage = "chi2cut-"+str(args.chi2cut[0])+"-"+str(args.chi2cut[1])+"_"+titlepage
            if args.mask is not None :
                titlepage = args.mask+"_"+titlepage
            model = plotframes.create_model(thespec, thezb)
            html_dir = os.path.join(webdir,"pix"+pixel)
            if not os.path.exists(html_dir) : 
                os.makedirs(html_dir)
                os.mkdir(html_dir+"/vignettes")
            
            plotframes.plotspectra(thespec, zcatalog=zbest, model_from_zcat=True, vidata=None, model=None, title=titlepage, html_dir=html_dir, is_coadded=True, mask_type=args.mask_type)
            for i_spec in range(thespec.num_spectra()) :
                saveplot = html_dir+"/vignettes/pix"+pixel+"_"+str(i_page)+"_"+str(i_spec)+".png"
                utils_specviewer.miniplot_spectrum(thespec, i_spec, model=model, saveplot=saveplot, smoothing = args.vignette_smoothing)
            nspec_done += thespec.num_spectra()
        
        # Stop running if needed, only once a full pixel is completed
        if args.nmax_spectra is not None :
            if nspec_done >= args.nmax_spectra :
                log.info(str(nspec_done)+" spectra done : no other pixel will be processed")
                break



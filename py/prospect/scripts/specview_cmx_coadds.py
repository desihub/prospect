# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
====================================
prospect.scripts.specview_cmx_coadds
====================================

Write static html files from "coadd" files in CMX data, sorted by exposure.
"""

import os, sys, glob
import argparse
import numpy as np
from astropy.table import Table, vstack
import astropy.io.fits

import desispec.io
from desiutil.log import get_logger
import desispec.spectra
import desispec.frame

from ..viewer import plotspectra
from ..myspecselect import myspecselect
from ..myspecupdate import myspecupdate
from ..utilities import specviewer_selection, match_redrock_zfit_to_spectra, match_zcat_to_spectra

# List of bad fibers in CMX data (see eg SB / KD emails 23-24/03/2020)
_bad_fibers_cmx = [
    [1000, 1250], # r2 amp A glowing spot CCD defect
    [2250, 2500], # r4 amp D serial register CTE problem
    [4500, 4750]  # b9 amp C readout salt-and-pepper noise
]

def _parse():

    parser = argparse.ArgumentParser(description='Create static html pages from CMX coadds, tile-based')
    parser.add_argument('--specprod_dir', help='Location of directory tree (data in specprod_dir/tiles/)', type=str)
    parser.add_argument('--cumulative', help='Use data from specprod_dir/tiles/cumulative/', action='store_true')
    parser.add_argument('--tile', help='Name of single tile to be processed',type=str, default=None)
    parser.add_argument('--tile_list', help='ASCII file providing list of tiles', type=str, default=None)
    parser.add_argument('--night', help='Filter night to be processed (night name can also be "deep")',type=str, default=None)
    parser.add_argument('--night_list', help='Filter set of nights to be processed (can include night name "deep")', type=str, default=None)
    parser.add_argument('--nspecperfile', help='Number of spectra in each html page', type=int, default=50)
    parser.add_argument('--webdir', help='Base directory for webpages', type=str)
    parser.add_argument('--nmax_spectra', help='Stop the production of HTML pages once a given number of spectra are done', type=int, default=None)
    parser.add_argument('--mask', help='Select only objects with a given target mask', type=str, default=None)
    parser.add_argument('--mask_type', help='Target mask type', type=str, default='CMX_TARGET')
    parser.add_argument('--snrcut', help='Select only objects in a given range for MEDIAN_CALIB_SNR_B+R+Z', nargs='+', type=float, default=None)
    parser.add_argument('--with_zcatalog', help='Include redshift fit results (zbest files)', action='store_true')
    parser.add_argument('--with_multiple_models', help='Display several models (requires full redrock outputs)', action='store_true')
    parser.add_argument('--petals', help='Select only a set of petals (labelled 0 to 9)', nargs='+', type=str, default=None)
    parser.add_argument('--clean_bad_fibers_cmx', help='Remove list of known bad fibers (CMX conditions)', action='store_true')
    parser.add_argument('--template_dir', help='Redrock template directory', type=str, default=None)
    parser.add_argument('--top_metadata', help='List of highlighted metadata', nargs='+', type=str, default=None)

    args = parser.parse_args()
    return args


def tile_db(specprod_dir, tile_subset=None, night_subset=None, petals=None, with_zcatalog=False, cumulative=False) :
    '''
    Returns [ {tile, night, petals} for all tile/night available in specprod_dir/tiles tree ],
        with b,r,z frames whose name matches frametype
        tile_subset : list; if None, all available tiles will be included in the list
        night_subset : list; if not None, only coadds from these nights are included
        petals (list of strings) : select only data from a set of petals
        with_zcatalog : select only petals with a zbest file
    '''

    if petals is None : petals = [str(i) for i in range(10)]
    tiles_db = list()
    basedir = os.path.join(specprod_dir,'tiles')
    if cumulative: basedir = os.path.join(basedir,'cumulative')
    for tile in os.listdir(basedir) :
        if tile_subset is not None and tile not in tile_subset : continue
        for night in os.listdir( os.path.join(basedir,tile) ) :
            if night_subset is not None and night not in night_subset : continue
            petals_avail = []
            for petal_num in petals :
                add_petal = True
                night_label = 'thru'+night if cumulative else night
                filename = "coadd-"+petal_num+"-"+tile+"-"+night_label+".fits"
                if not os.path.isfile( os.path.join(basedir,tile,night,filename) ) : add_petal=False
                if with_zcatalog :
                    filename = "zbest-"+petal_num+"-"+tile+"-"+night_label+".fits"
                    if not os.path.isfile( os.path.join(basedir,tile,night,filename) ) : add_petal=False
                if add_petal :
                    petals_avail.append(petal_num)
            if len(petals_avail) > 0 :
                tiles_db.append( { 'tile':tile, 'night':night, 'petals':petals_avail} )

    return tiles_db


def page_subset_tile(fdir, tile_db_subset, html_dir, titlepage_prefix, mask, log, nspecperfile, snr_cut, with_zcatalog=False, template_dir=None, clean_bad_fibers_cmx=False, with_multiple_models=False, mask_type='CMX_TARGET', top_metadata=None) :
    '''
    Running prospect from coadds.
    '''

    tile = tile_db_subset['tile']
    night = tile_db_subset['night']
    night_label = 'thru'+night if 'cumulative' in fdir else night
    nspec_done = 0
    all_spectra = None
    log.info("Tile "+tile+": reading coadds from night "+night)
    if with_multiple_models :
        rrtables = []

    for petal_num in tile_db_subset['petals'] :
        fname = os.path.join(fdir,"coadd-"+petal_num+"-"+tile+'-'+night_label+".fits")
        spectra = desispec.io.read_spectra(fname)
        # Filtering
        if (mask != None) or (snr_cut != None) :
            spectra = specviewer_selection(spectra, log=log,
                        mask=mask, mask_type=mask_type, snr_cut=snr_cut, with_dirty_mask_merge=True)
            if spectra == 0 : continue
        # Display multiple models: requires redrock catalog
        # Hardcoded: display up to 4th best fit (=> need 5 best fits in redrock table)
        if with_multiple_models :
            fname = os.path.join(fdir,'redrock-'+petal_num+'-'+tile+'-'+night_label+".h5")
            rr_table = match_redrock_zfit_to_spectra(fname, spectra)
            rrtables.append(rr_table)
        # Merge
        if all_spectra is None :
            all_spectra = spectra
        else :
            # NB update() does not copy scores. Score-based filtering (SNR) was done before.
            all_spectra = myspecupdate(all_spectra, spectra)

    if all_spectra is None :
        log.info("Tile "+tile+" - night "+night+": no spectra.")
        return 0
    else :
        clean_fiberstatus=True if 'FIBERSTATUS' in all_spectra.fibermap.keys() else False
        fibers = None
        if clean_bad_fibers_cmx :
            fibers = np.arange(5000)
            for cut_fiber_range in _bad_fibers_cmx :
                fibers = fibers[ ( (fibers < cut_fiber_range[0]) | (fibers > cut_fiber_range[1]) )]
        all_spectra = myspecselect(all_spectra,
                            clean_fiberstatus=clean_fiberstatus, fibers=fibers, remove_scores=True)
        if all_spectra is None : return 0

    # zcatalog (adding redrock version)
    if with_zcatalog :
        ztables = []
        for petal_num in tile_db_subset['petals'] :
            fname = os.path.join(fdir,"zbest-"+petal_num+"-"+tile+'-'+night_label+".fits")
            the_ztable = Table.read(fname,'ZBEST')
            hdulist = astropy.io.fits.open(fname)
            the_ztable['RRVER'] = hdulist[hdulist.index_of('PRIMARY')].header['RRVER']
            ztables.append(the_ztable)
        zcat = vstack(ztables)
    else : zcat = None

    if with_multiple_models :
        rrtable = vstack(rrtables)
        num_approx_fits = 4 # TODO settle option/unhardcode
        with_full_2ndfit = True # TODO settle option/unhardcode
    else :
        rrtable = None
        num_approx_fits = None
        with_full_2ndfit = False

    # Create several html pages : sort by targetid
    nspec_tile = all_spectra.num_spectra()
    log.info("Tile "+tile+" - night "+night+": "+str(nspec_tile)+" exposure-coadded spectra")
    sort_indices = np.argsort(all_spectra.fibermap["TARGETID"])
    nbpages = int(np.ceil((nspec_tile/nspecperfile)))
    for i_page in range(1,1+nbpages) :

        log.info(" * Page "+str(i_page)+" / "+str(nbpages))
        the_indices = sort_indices[(i_page-1)*nspecperfile:i_page*nspecperfile]
        thespec = myspecselect(all_spectra, indices=the_indices, remove_scores=True)
        the_zcat, kk = match_zcat_to_spectra(zcat, thespec)
        if with_multiple_models :
            the_rrtable, kk = match_zcat_to_spectra(rrtable, thespec)
        else :
            the_rrtable = None

        titlepage = titlepage_prefix+"_"+str(i_page)
        plotspectra(thespec, with_noise=True, zcatalog=the_zcat,
                    title=titlepage, html_dir=html_dir, mask_type=mask_type, with_thumb_only_page=True,
                    template_dir=template_dir, redrock_cat=the_rrtable, num_approx_fits=num_approx_fits,
                    with_full_2ndfit=with_full_2ndfit, top_metadata=top_metadata)
    nspec_done += nspec_tile

    return nspec_done



def main():
    args = _parse()
    log = get_logger()
    webdir = args.webdir
    if ( [args.tile, args.tile_list] ).count(None) != 1 :
        log.info("Specview_cmx_coadds: Wrong set of input tiles. Exiting")
        return 0
    if args.night is not None and args.night_list is not None :
        log.info("Specview_cmx_coadds: Wrong set of input nights. Exiting")
        return 0
    
    # Logistics : list of "subsets" to process
    if args.tile_list is not None :
        tile_subset = np.loadtxt(args.tile_list, dtype=str, comments='#')
    if args.tile is not None :
        tile_subset = [ args.tile ]
    if args.night is not None :
        night_subset = [ args.night ]
    elif args.night_list is not None :
        night_subset = np.loadtxt(args.night_list, dtype=str, comments='#')
    else : night_subset = None
    subset_db = tile_db(args.specprod_dir, tile_subset=tile_subset, night_subset=night_subset,
                        petals=args.petals, with_zcatalog=args.with_zcatalog, cumulative=args.cumulative)
    tmplist = [ x['tile'] for x in subset_db ]
    missing_tiles = [ x for x in tile_subset if x not in tmplist ]
    for x in missing_tiles : log.info("Missing tile, cannot be processed: "+x)
    log.info(str(len(subset_db))+" tile/night subsets to be processed")

    # Main loop on subsets
    nspec_done = 0
    for the_subset in subset_db :

        log.info("Working on tile "+the_subset['tile']+" - night "+the_subset['night'])
        if args.cumulative:
            fdir = os.path.join( args.specprod_dir, 'tiles', 'cumulative', the_subset['tile'], the_subset['night'] )
        else:
            fdir = os.path.join( args.specprod_dir, 'tiles', the_subset['tile'], the_subset['night'] )
        html_dir = os.path.join(webdir,"tiles",the_subset['tile'],the_subset['night'])
        titlepage_prefix = "tile"+the_subset['tile']+"_night"+the_subset['night']
        if args.mask != None :
            html_dir = os.path.join(webdir, "tiles_"+args.mask, the_subset['tile'],the_subset['night'])
            titlepage_prefix = args.mask+"_"+titlepage_prefix
        if not os.path.exists(html_dir) :
            os.makedirs(html_dir)

        nspec_added = page_subset_tile(fdir, the_subset, html_dir, titlepage_prefix, args.mask, log, args.nspecperfile, args.snrcut, with_zcatalog=args.with_zcatalog, template_dir=args.template_dir, clean_bad_fibers_cmx=args.clean_bad_fibers_cmx, with_multiple_models=args.with_multiple_models, mask_type=args.mask_type, top_metadata=args.top_metadata)

        # Stop running if needed, only once a full exposure is completed
        nspec_done += nspec_added
        if args.nmax_spectra is not None :
            if nspec_done >= args.nmax_spectra :
                log.info(str(nspec_done)+" spectra done: no other exposure will be processed")
                break

    log.info("End of specview_cmx_coadds script.")
    return 0

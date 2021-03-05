# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
=====================================
prospect.scripts.specview_cmx_targets
=====================================

Write static html files from a set of targets, using tile-based coadds in CMX data.
"""

import os
import argparse
import numpy as np

import desispec.io
from desiutil.log import get_logger

from ..viewer import plotspectra
from ..myspecselect import myspecselect
from ..utilities import make_targetdict, load_spectra_zcat_from_targets  #, match_zcat_to_spectra

def _parse():

    parser = argparse.ArgumentParser(description='Create static html pages from a set of targets, using CMX tile-based coadds')
    parser.add_argument('--specprod_dir', help='Location of directory tree (data in specprod_dir/tiles/)', type=str)
    parser.add_argument('--target_list', help='ASCII file providing the list of targetids', type=str)
    parser.add_argument('--tiles', help='Name of tile[s] to be processed (avoids to scan all tiles)', nargs='+', type=str, default=None)
    parser.add_argument('--nights', help='Name of night[s] to be processed (avoids to scan all nights)', nargs='+', type=str, default=None)
    parser.add_argument('--nspecperfile', help='Number of spectra in each html page', type=int, default=50)
    parser.add_argument('--webdir', help='Base directory for webpages', type=str)
    parser.add_argument('--titlepage_prefix', help='Prefix for webpage title', type=str, default='targetlist')
    parser.add_argument('--with_multiple_models', help='Display several models (requires full redrock outputs)', action='store_true')
    parser.add_argument('--template_dir', help='Redrock template directory', type=str, default=None)
    parser.add_argument('--countdown', help='Insert a countdown widget (in minutes)', type=int, default=-1)

    args = parser.parse_args()
    return args


def main():
    args = _parse()
    log = get_logger()

    tile_dir = os.path.join(args.specprod_dir,'tiles')
    if args.tiles :
        tile_list = args.tiles
    else :
        tile_list = os.listdir(tile_dir)
    obs_db = make_targetdict(tile_dir, tiles=tile_list, nights=args.nights)

    targetids = np.loadtxt(args.target_list, dtype='int64', comments='#')
    log.info(str(len(targetids))+" targets provided.")

    if args.with_multiple_models :
        spectra, zcat, rrtable = load_spectra_zcat_from_targets(targetids, tile_dir, obs_db, with_redrock=True)
    else :
        spectra, zcat = load_spectra_zcat_from_targets(targetids, tile_dir, obs_db, with_redrock=False)

    # TODO? this may be put in a standalone fct, avoid code duplicate with other script
    # Create several html pages : sort by targetid
    nspec = spectra.num_spectra()
    log.info(str(nspec)+" spectra obtained.")
    sort_indices = np.argsort(spectra.fibermap["TARGETID"])
    nbpages = int(np.ceil((nspec/args.nspecperfile)))
    for i_page in range(1,1+nbpages) :

        log.info(" * Page "+str(i_page)+" / "+str(nbpages))
        the_indices = sort_indices[(i_page-1)*args.nspecperfile:i_page*args.nspecperfile]
        thespec, kept_ind = myspecselect(spectra, indices=the_indices, remove_scores=True, output_indices=True)
        the_zcat = zcat[kept_ind]
        if not np.array_equal(the_zcat['TARGETID'], thespec.fibermap['TARGETID']) :
            raise RuntimeError("targetids do not match between spec and zcat")
        #the_zcat, kk = match_zcat_to_spectra(zcat, thespec)
        if args.with_multiple_models :
            the_rrtable = rrtable[kept_ind]
            #the_rrtable, kk = match_zcat_to_spectra(rrtable, thespec)
            num_approx_fits = 4 # TODO settle option
            with_full_2ndfit = True # TODO settle option
        else :
            the_rrtable = None
            num_approx_fits = None
            with_full_2ndfit = False

        titlepage = args.titlepage_prefix+"_"+str(i_page)
        plotspectra(thespec, with_noise=True, zcatalog=the_zcat,
                    title=titlepage, html_dir=args.webdir, mask_type='CMX_TARGET', with_thumb_only_page=True,
                    template_dir=args.template_dir, redrock_cat=the_rrtable, num_approx_fits=num_approx_fits,
                    with_full_2ndfit=with_full_2ndfit, vi_countdown=args.countdown)

    log.info("End of specview_cmx_targets script.")
    return 0

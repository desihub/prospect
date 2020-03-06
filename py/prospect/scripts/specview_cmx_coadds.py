"""
prospect.scripts.specview_cmx_frames
===================================

Write static html files from "coadd" files in CMX data
"""

import os, sys, glob
import argparse
import numpy as np
from astropy.table import Table, vstack

import desispec.io
from desiutil.log import get_logger
from desitarget.targetmask import desi_mask
from desitarget.cmx.cmx_targetmask import cmx_mask
import desispec.spectra
import desispec.frame

from prospect import plotframes
from prospect import utils_specviewer
from prospect import myspecselect


def parse() :

    parser = argparse.ArgumentParser(description='Create static html pages from CMX coadds, tile-based')
    parser.add_argument('--specprod_dir', help='Location of directory tree (data in specprod_dir/tiles/)', type=str)
    parser.add_argument('--tile', help='Name of single tile to be processed',type=str, default=None)
    parser.add_argument('--tile_list', help='ASCII file providing list of tiles', type=str, default=None)
    parser.add_argument('--nspecperfile', help='Number of spectra in each html page', type=int, default=50)
    parser.add_argument('--webdir', help='Base directory for webpages', type=str)
    parser.add_argument('--nmax_spectra', help='Stop the production of HTML pages once a given number of spectra are done', type=int, default=None)
    parser.add_argument('--mask', help='Select only objects with a given CMX_TARGET target mask', type=str, default=None)
    parser.add_argument('--snrcut', help='Select only objects in a given range for MEDIAN_CALIB_SNR_B+R+Z', nargs='+', type=float, default=None)
    parser.add_argument('--with_zcatalog', help='Include redshift fit results (zbest files)', action='store_true')
    parser.add_argument('--petals', help='Select only a set of petals (labelled 0 to 9)', nargs='+', type=str, default=None)

    args = parser.parse_args()
    return args


def tile_db(specprod_dir, tile_subset=None, night_subset=None, petals=None, with_zcatalog=False) :
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
    for tile in os.listdir( os.path.join(specprod_dir,'tiles') ) :
        if tile_subset is not None and tile not in tile_subset : continue
        for night in os.listdir( os.path.join(specprod_dir,'tiles',tile) ) :
            if night_subset is not None and night not in night_subset : continue
            petals_avail = []
            for petal_num in petals :
                add_petal = True
                filename = "coadd-"+petal_num+"-"+tile+"-"+night+".fits"
                if not os.path.isfile( os.path.join(specprod_dir,'tiles',tile,night,filename) ) : add_petal=False
                if with_zcatalog :
                    filename = "zbest-"+petal_num+"-"+tile+"-"+night+".fits"
                    if not os.path.isfile( os.path.join(specprod_dir,'tiles',tile,night,filename) ) : add_petal=False
                if add_petal :
                    petals_avail.append(petal_num)
            if len(petals_avail) > 0 :
                tiles_db.append( { 'tile':tile, 'night':night, 'petals':petals_avail} )
                    
    return tiles_db
    

def page_subset_tile(fdir, tile_db_subset, html_dir, titlepage_prefix, mask, log, nspecperfile, snr_cut, with_zcatalog=False) :
    '''
    Running prospect from coadds.
    '''
    
    tile = tile_db_subset['tile']
    night = tile_db_subset['night']
    nspec_done = 0
    all_spectra = None
    log.info("Tile "+tile+" : reading coadds from night "+night)
    for petal_num in tile_db_subset['petals'] :
        fname = os.path.join(fdir,"coadd-"+petal_num+"-"+tile+'-'+night+".fits")
        spectra = desispec.io.read_spectra(fname)
        ### TMP TRICK (?) : need to have (DUMMY) exposures in fibermap, otherwise spectra.update() crashes !
        if not('EXPID' in spectra.fibermap.keys()) :
                spectra.fibermap['EXPID'] = spectra.fibermap['FIBER']*0
        ### END TMP TRICK
        # Filtering
        if (mask != None) or (snr_cut != None) :
            spectra = utils_specviewer.specviewer_selection(spectra, log=log,
                        mask=mask, mask_type='CMX_TARGET', snr_cut=snr_cut)
            if spectra == 0 : continue
        # Merge
        if all_spectra is None :
            all_spectra = spectra
        else :
            all_spectra.update(spectra) # NB update() does not copy scores. Filtering was done before.
                
    if all_spectra is None : 
        log.info("Tile "+tile+" - night "+night+": no spectra.")
        return 0

    # zcatalog
    if with_zcatalog :
        ztables = []
        for petal_num in tile_db_subset['petals'] :
            fname = os.path.join(fdir,"zbest-"+petal_num+"-"+tile+'-'+night+".fits")
            ztables.append(Table.read(fname,'ZBEST'))
        zcat = vstack(ztables)
    else : zcat = None
    
    # Create several html pages : sort by targetid
    nspec_tile = all_spectra.num_spectra()
    log.info("Tile "+tile+" - night "+night+": "+str(nspec_tile)+" exposure-coadded spectra")
    sort_indices = np.argsort(all_spectra.fibermap["TARGETID"])
    nbpages = int(np.ceil((nspec_tile/nspecperfile)))
    for i_page in range(1,1+nbpages) :

        log.info(" * Page "+str(i_page)+" / "+str(nbpages))
        the_indices = sort_indices[(i_page-1)*nspecperfile:i_page*nspecperfile]            
        thespec = myspecselect.myspecselect(all_spectra, indices=the_indices)
        titlepage = titlepage_prefix+"_"+str(i_page)
        plotframes.plotspectra(thespec, with_noise=True, with_coaddcam=False, is_coadded=True, zcatalog=zcat,
                    title=titlepage, html_dir=html_dir, mask_type='CMX_TARGET', with_thumb_only_page=True)
    nspec_done += nspec_tile
        
    return nspec_done



def main(args) :
    
    log = get_logger()
    webdir = args.webdir
    if ( [args.tile, args.tile_list] ).count(None) != 1 :
        log.info("Specview_cmx_coadds : Wrong set of input tiles. Exiting")
        return 0
    
    # Logistics : list of "subsets" to process        
    if args.tile_list is not None :
        tile_subset = np.loadtxt(args.tile_list, dtype=str, comments='#')
    if args.tile is not None :
        tile_subset = [ args.tile ]
    subset_db = tile_db(args.specprod_dir, tile_subset=tile_subset, petals=args.petals, with_zcatalog=args.with_zcatalog)
    tmplist = [ x['tile'] for x in subset_db ]
    missing_tiles = [ x for x in tile_subset if x not in tmplist ]
    for x in missing_tiles : log.info("Missing tile, cannot be processed : "+x)
    log.info(str(len(subset_db))+" tiles [exposures] to be processed")
    
    # Main loop on subsets
    nspec_done = 0
    for the_subset in subset_db :
        
        log.info("Working on tile "+the_subset['tile']+" - night "+the_subset['night'])
        fdir = os.path.join( args.specprod_dir, 'tiles', the_subset['tile'], the_subset['night'] )
        html_dir = os.path.join(webdir,"tiles",the_subset['tile'],the_subset['night'])
        titlepage_prefix = "tile"+the_subset['tile']+"_night"+the_subset['night']
        if args.mask != None :
            html_dir = os.path.join(webdir, "tiles_"+args.mask, the_subset['tile'],the_subset['night'])
            titlepage_prefix = args.mask+"_"+titlepage_prefix
        if not os.path.exists(html_dir) : 
            os.makedirs(html_dir)
        
        nspec_added = page_subset_tile(fdir, the_subset, html_dir, titlepage_prefix, args.mask, log, args.nspecperfile, args.snrcut, with_zcatalog=args.with_zcatalog)
                    
        # Stop running if needed, only once a full exposure is completed
        nspec_done += nspec_added
        if args.nmax_spectra is not None :
            if nspec_done >= args.nmax_spectra :
                log.info(str(nspec_done)+" spectra done : no other exposure will be processed")
                break

        return 0

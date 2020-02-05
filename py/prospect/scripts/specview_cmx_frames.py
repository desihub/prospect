"""
prospect.scripts.specview_cmx_frames
===================================

Write static html files from "[s/c/..]frame" files in CMX data
Don't know if this will be used for SV (once frames are superseeded by spectra)
"""

import os, sys, glob
import argparse
import numpy as np

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

    parser = argparse.ArgumentParser(description='Create static html pages from CMX frames, tile- or exposure-based')
    parser.add_argument('--specprod_dir', help='Location of directory tree (data in specprod_dir/exposures/)', type=str)
    parser.add_argument('--exposure', help='Name of single exposure to be processed',type=str, default=None)
    parser.add_argument('--exposure_list', help='ASCII file providing list of exposures', type=str, default=None)
    parser.add_argument('--tile', help='Name of single tile to be processed',type=str, default=None)
    parser.add_argument('--tile_list', help='ASCII file providing list of tiles', type=str, default=None)
    parser.add_argument('--nspecperfile', help='Number of spectra in each html page', type=int, default=50)
    parser.add_argument('--webdir', help='Base directory for webpages', type=str)
    parser.add_argument('--nmax_spectra', help='Stop the production of HTML pages once a given number of spectra are done', type=int, default=None)
    parser.add_argument('--frametype', help='Input frame category (currently sframe/cframe supported)', type=str, default='cframe')
    parser.add_argument('--mask', help='Select only objects with a given CMX_TARGET target mask', type=str, default=None)
    args = parser.parse_args()
    return args


def exposure_db(specprod_dir, frametype='cframe', expo_subset=None) :
    '''
    Returns list of {exposure, night, spectrographs} available in specprod_dir/exposures tree, 
        with b,r,z frames whose name matches frametype
        spectros = list of spectrograph numbers from 0 to 9 (only those with 3 bands available are kept)
        expo_subset : list; if None, all available exposures will be included in the list
    '''
    expo_db = list()
    if expo_subset is not None :
        # List of exposures (eg from an ASCII file) does not need zero-padding
        expo_subset = [ x.rjust(8, "0") for x in expo_subset]
    for night in os.listdir( os.path.join(specprod_dir,'exposures') ) :
        for expo in os.listdir( os.path.join(specprod_dir,'exposures',night) ) :
            if not os.path.isdir( os.path.join(specprod_dir,'exposures',night,expo) ) : continue
            if expo_subset is not None :
                if expo not in expo_subset : continue
            files_avail = os.listdir( os.path.join(specprod_dir,'exposures',night,expo) )
            spectro_avail = []
            for spectro_num in [str(i) for i in range(10)] :
                files_check = [ frametype+"-"+band+spectro_num+"-"+expo+".fits" for band in ['b','r','z'] ]
                if all(x in files_avail for x in files_check) :
                    spectro_avail.append(spectro_num)
            if len(spectro_avail) > 0 :
                expo_db.append( {'exposure':expo, 'night':night, 'spectrographs':spectro_avail} )
    return expo_db


def tile_db(specprod_dir, frametype='cframe', tile_subset=None, night_subset=None) :
    '''
    Returns list of {tile, night, expo, spectros} available in specprod_dir/tiles tree, 
        with b,r,z frames whose name matches frametype
        spectros = list of spectrograph numbers from 0 to 9 (only those with 3 bands available are kept)
        tile_subset : list; if None, all available tiles will be included in the list
        night_subset : list; if not None, only frames from these nights are included
    '''
    tiles_db = list()
    for tile in os.listdir( os.path.join(specprod_dir,'tiles') ) :
        if tile not in tile_subset : 
            continue
        for night in os.listdir( os.path.join(specprod_dir,'tiles',tile) ) :
            if night_subset is not None :
                if night not in night_subset : continue
            files_avail = os.listdir( os.path.join(specprod_dir,'tiles',tile,night) )
            files_avail = [ x for x in files_avail if frametype in x and x[-5:]==".fits"]
            expos_avail = list( set([ x[-13:-5] for x in files_avail ]) ) # This assumes files are of the type *-expid.fits

            for expo in expos_avail :
                spectro_avail = []
                for spectro_num in [str(i) for i in range(10)] :
                    files_check = [ frametype+"-"+band+spectro_num+"-"+expo+".fits" for band in ['b','r','z'] ]
                    if all(x in files_avail for x in files_check) :
                        spectro_avail.append(spectro_num)
                if len(spectro_avail) > 0 :
                    tiles_db.append( { 'tile':tile, 'night':night, 'exposure':expo, 'spectrographs':spectro_avail} )
    return tiles_db


def page_subset(fdir, exposure, frametype, spectrographs, html_dir, titlepage_prefix, mask, log, nspecperfile) :
    '''
    Running prospect from frames : loop over spectrographs for a given exposure
    '''
    
    nspec_done = 0
    for spectrograph_num in spectrographs :
        frames = [ desispec.io.read_frame(os.path.join(fdir,frametype+"-"+band+spectrograph_num+"-"+exposure+".fits")) for band in ['b','r','z'] ]
        spectra = utils_specviewer.frames2spectra(frames)
        # Selection
        if mask != None :
            spectra = utils_specviewer.specviewer_selection(spectra, log=log,
                        mask=mask, mask_type='CMX_TARGET')
            if spectra == 0 : continue

        # Handle several html pages per exposure - sort by fiberid
        nspec_expo = spectra.num_spectra()
        log.info("Spectrograph number "+spectrograph_num+" : "+str(nspec_expo)+" spectra")
        sort_indices = np.argsort(frames[0].fibermap["FIBER"])
        nbpages = int(np.ceil((nspec_expo/nspecperfile)))
        for i_page in range(1,1+nbpages) :

            log.info(" * Page "+str(i_page)+" / "+str(nbpages))
            the_indices = sort_indices[(i_page-1)*nspecperfile:i_page*nspecperfile]            
            thespec = myspecselect.myspecselect(spectra, indices=the_indices)
            titlepage = titlepage_prefix+"_spectro"+spectrograph_num+"_"+str(i_page)
            plotframes.plotspectra(thespec, with_noise=True, with_coaddcam=True, is_coadded=False, 
                        title=titlepage, html_dir=html_dir, mask_type='CMX_TARGET', with_thumb_only_page=True)
        nspec_done += nspec_expo
        
    return nspec_done


def main(args) :
    
    log = get_logger()
    webdir = args.webdir
    assert args.frametype in ['sframe', 'cframe']
    if ( [args.exposure_list, args.exposure, args.tile, args.tile_list] ).count(None) != 3 :
        log.info("Specview_cmx_frames : Wrong set of input tiles/exposures. Exiting")
        return 0
    
    page_sorting = "tile"
    if (args.exposure_list is not None) or (args.exposure is not None) : page_sorting = "exposure"

    # Logistics : list of "subsets" to process
    if page_sorting == "exposure" :
        expo_subset = None
        if args.exposure_list is not None :
            expo_subset = np.loadtxt(args.exposure_list, dtype=str, comments='#')
        if args.exposure is not None :
            expo_subset = [ args.exposure ]
        subset_db = exposure_db(args.specprod_dir, frametype=args.frametype, expo_subset=expo_subset)

        if expo_subset is not None :
            tmplist = [ x['exposure'] for x in subset_db ]
            missing_expos = [ x for x in expo_subset if x.rjust(8, "0") not in tmplist ]
            for x in missing_expos : log.info("Missing exposure, cannot be processed : "+x)
        log.info(str(len(subset_db))+" exposures to be processed")    
    if page_sorting == "tile" :
        tile_subset = None
        if args.tile_list is not None :
            tile_subset = np.loadtxt(args.tile_list, dtype=str, comments='#')
        if args.tile is not None :
            tile_subset = [ args.tile ]
        subset_db = tile_db(args.specprod_dir, frametype=args.frametype, tile_subset=tile_subset) # TODO select night subset

        if tile_subset is not None :
            tmplist = [ x['tile'] for x in subset_db ]
            missing_tiles = [ x for x in tile_subset if x not in tmplist ]
            for x in missing_tiles : log.info("Missing tile, cannot be processed : "+x)
        log.info(str(len(subset_db))+" (tiles/exposures) to be processed")

    # Main loop on subsets
    nspec_done = 0
    for the_subset in subset_db :
        
        if page_sorting == 'exposure' :
            log.info("Working on exposure "+the_subset['exposure'])
            fdir = os.path.join( args.specprod_dir, 'exposures', the_subset['night'], the_subset['exposure'] )
            html_dir = os.path.join(webdir,"exposures_all",the_subset['exposure'])
            titlepage_prefix = "expo"+the_subset['exposure']
            if args.mask != None :
                html_dir = os.path.join(webdir, "exposures_"+args.mask, the_subset['exposure'])
                titlepage_prefix = args.mask+"_"+titlepage_prefix
        if page_sorting == 'tile' :
            log.info("Working on tile "+the_subset['tile']+" (exposure "+the_subset['exposure']+")")
            fdir = os.path.join( args.specprod_dir, 'tiles', the_subset['tile'], the_subset['night'] )
            html_dir = os.path.join(webdir,"tiles_all",the_subset['tile'])
            titlepage_prefix = "tile"+the_subset['tile']+"_expo"+the_subset['exposure']
            if args.mask != None :
                html_dir = os.path.join(webdir, "tiles_"+args.mask, the_subset['tile'])
                titlepage_prefix = args.mask+"_"+titlepage_prefix        
        if not os.path.exists(html_dir) : 
            os.makedirs(html_dir)
        
        nspec_added = page_subset(fdir, the_subset['exposure'], args.frametype, the_subset['spectrographs'], html_dir, titlepage_prefix, args.mask, log, args.nspecperfile)
                
        # Stop running if needed, only once a full exposure is completed
        nspec_done += nspec_added
        if args.nmax_spectra is not None :
            if nspec_done >= args.nmax_spectra :
                log.info(str(nspec_done)+" spectra done : no other exposure will be processed")
                break

        return 0



"""
prospect.scripts.specview_cmx_frames
===================================

Write static html files from "[s/c/..]frame" files in CMX data
Don't know if this will be used for SV (once frames are superseeded by spectra)
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

    parser = argparse.ArgumentParser(description='Create static html pages from CMX frames, tile- or exposure-based')
    parser.add_argument('--specprod_dir', help='Location of directory tree (data in specprod_dir/exposures/)', type=str)
    parser.add_argument('--exposure', help='Name of single exposure to be processed',type=str, default=None)
    parser.add_argument('--exposure_list', help='ASCII file providing list of exposures', type=str, default=None)
    parser.add_argument('--tile', help='Name of single tile to be processed',type=str, default=None)
    parser.add_argument('--tile_list', help='ASCII file providing list of tiles', type=str, default=None)
    parser.add_argument('--sort_exposures', help='If tile-based, will still sort pages by exposures/petals', action='store_true')
    parser.add_argument('--nspecperfile', help='Number of spectra in each html page', type=int, default=50)
    parser.add_argument('--webdir', help='Base directory for webpages', type=str)
    parser.add_argument('--nmax_spectra', help='Stop the production of HTML pages once a given number of spectra are done', type=int, default=None)
    parser.add_argument('--frametype', help='Input frame category (currently sframe/cframe supported)', type=str, default='cframe')
    parser.add_argument('--mask', help='Select only objects with a given CMX_TARGET target mask', type=str, default=None)
    parser.add_argument('--snrcut', help='Select only objects in a given range for MEDIAN_CALIB_SNR_B+R+Z', nargs='+', type=float, default=None)
    parser.add_argument('--with_zcatalog', help='If tile-based, will include redshift fit results (zbest files)', action='store_true')
    parser.add_argument('--petals', help='Select only a set of petals (labelled 0 to 9)', nargs='+', type=str, default=None)

    args = parser.parse_args()
    return args


def exposure_db(specprod_dir, frametype='cframe', expo_subset=None, petals=None) :
    '''
    Returns list of {exposure, night, petals} available in specprod_dir/exposures tree, 
        with b,r,z frames whose name matches frametype
        petals = list of petal numbers from 0 to 9 (only those with 3 bands available are kept)
        expo_subset : list; if None, all available exposures will be included in the list
        petals (list of strings) : select only data from a set of petals
    '''
    if petals is None : petals = [str(i) for i in range(10)]
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
            petals_avail = []
            for petal_num in petals :
                files_check = [ frametype+"-"+band+petal_num+"-"+expo+".fits" for band in ['b','r','z'] ]
                if all(x in files_avail for x in files_check) :
                    petals_avail.append(petal_num)
            if len(petals_avail) > 0 :
                expo_db.append( {'exposure':expo, 'night':night, 'petals':petals_avail} )
    return expo_db


def tile_db(specprod_dir, frametype='cframe', tile_subset=None, night_subset=None, merge_exposures=False, petals=None) :
    '''
    Returns [ {tile, night, expo, petals} for all tile/expos available in specprod_dir/tiles tree ], 
        with b,r,z frames whose name matches frametype
        tile_subset : list; if None, all available tiles will be included in the list
        night_subset : list; if not None, only frames from these nights are included
        merge_exposures : if True, returns [ [ {tile, night, expo, petals} for all expos ] for all tiles ]
        petals (list of strings) : select only data from a set of petals, from 0 to 9 (only those with 3 band data available are kept)
    '''
    if petals is None : petals = [str(i) for i in range(10)]
    tiles_db = list()
    for tile in os.listdir( os.path.join(specprod_dir,'tiles') ) :
        if tile_subset is not None and if tile not in tile_subset : 
            continue
        if merge_exposures : tile_subdb = { 'tile':tile, 'db_subset':list() }
        for night in os.listdir( os.path.join(specprod_dir,'tiles',tile) ) :
            if night_subset is not None :
                if night not in night_subset : continue
            files_avail = os.listdir( os.path.join(specprod_dir,'tiles',tile,night) )
            files_avail = [ x for x in files_avail if frametype in x and x[-5:]==".fits"]
            expos_avail = list( set([ x[-13:-5] for x in files_avail ]) ) # This assumes files are of the type *-expid.fits

            for expo in expos_avail :
                petals_avail = []
                for petal_num in petals :
                    files_check = [ frametype+"-"+band+petal_num+"-"+expo+".fits" for band in ['b','r','z'] ]
                    if all(x in files_avail for x in files_check) :
                        petals_avail.append(petal_num)
                if len(petals_avail) > 0 :
                    if merge_exposures :
                        tile_subdb['db_subset'].append( { 'tile':tile, 'night':night, 'exposure':expo, 'petals':petals_avail} )
                    else :
                        tiles_db.append( { 'tile':tile, 'night':night, 'exposure':expo, 'petals':petals_avail} )
        if merge_exposures and ( len(tile_subdb['db_subset'])>0 ) : tiles_db.append(tile_subdb)
    return tiles_db    
    

def page_subset_expo(fdir, exposure, frametype, petals, html_dir, titlepage_prefix, mask, log, nspecperfile, snr_cut) :
    '''
    Running prospect from frames : loop over petals for a given exposure
    '''
    
    nspec_done = 0
    for petal_num in petals :
        frames = [ desispec.io.read_frame(os.path.join(fdir,frametype+"-"+band+petal_num+"-"+exposure+".fits")) for band in ['b','r','z'] ]
        spectra = utils_specviewer.frames2spectra(frames, with_scores=True)
        # Selection
        if (mask != None) or (snr_cut != None) :
            spectra = utils_specviewer.specviewer_selection(spectra, log=log,
                        mask=mask, mask_type='CMX_TARGET', snr_cut=snr_cut)
            if spectra == 0 : continue

        # Handle several html pages per exposure - sort by fiberid
        nspec_expo = spectra.num_spectra()
        log.info("Petal number "+petal_num+" : "+str(nspec_expo)+" spectra")
        sort_indices = np.argsort(spectra.fibermap["FIBER"])
        nbpages = int(np.ceil((nspec_expo/nspecperfile)))
        for i_page in range(1,1+nbpages) :

            log.info(" * Page "+str(i_page)+" / "+str(nbpages))
            the_indices = sort_indices[(i_page-1)*nspecperfile:i_page*nspecperfile]            
            thespec = myspecselect.myspecselect(spectra, indices=the_indices)
            titlepage = titlepage_prefix+"_petal"+petal_num+"_"+str(i_page)
            plotframes.plotspectra(thespec, with_noise=True, with_coaddcam=True, is_coadded=False, 
                        title=titlepage, html_dir=html_dir, mask_type='CMX_TARGET', with_thumb_only_page=True)
        nspec_done += nspec_expo
        
    return nspec_done

def page_subset_tile(fdir, tile_db_subset, frametype, html_dir, titlepage_prefix, mask, log, nspecperfile, snr_cut, with_zcatalog=False) :
    '''
    Running prospect from frames : tile-based, do not separate pages per exposure.
        tile_db_subset : subset of tile_db, all with the same tile
    '''
    
    tile = tile_db_subset['tile']
    nspec_done = 0
    all_spectra = None
    for the_subset in tile_db_subset['db_subset'] :
        assert (the_subset['tile']==tile)
        log.info("Tile "+tile+" : reading frames from exposure "+the_subset['exposure'])
        for petal_num in the_subset['petals'] :
            frames = [ desispec.io.read_frame(os.path.join(fdir, the_subset['night'],frametype+"-" + band + petal_num + "-" + the_subset['exposure'] + ".fits")) for band in ['b','r','z'] ]
            ### TMP TRICK (?) : need to have exposures in fibermaps, otherwise spectra.update() crashes !
            for fr in frames :
                if not('EXPID' in fr.fibermap.keys()) :
                    fr.fibermap['EXPID'] = fr.fibermap['FIBER']
                    for i in range(len(fr.fibermap)) : fr.fibermap['EXPID'][i] = the_subset['exposure']
            ### END TMP TRICK
            ### OTHER TRICK : need resolution data in spectra to pass coadd fct (could be changed...)
            spectra = utils_specviewer.frames2spectra(frames, with_scores=True, with_resolution_data=True)
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
        log.info("Tile "+tile+" : no spectra !")
        return 0
    # Exposure-coadd
    all_spectra = utils_specviewer.coadd_targets(all_spectra)
    # zcatalog
    if with_zcatalog :
        zcat_files = glob.glob(fdir+"/"+the_subset['night']+"/zbest*.fits") # Probably TMP ... ?
        ztables = []
        for f in zcat_files : 
            ztables.append(Table.read(f,'ZBEST'))
        zcat = vstack(ztables)
    else : zcat = None
    # Handle several html pages per exposure - sort by targetid
    nspec_tile = all_spectra.num_spectra()
    log.info("Tile "+tile+" : "+str(nspec_tile)+" exposure-coadded spectra")
    sort_indices = np.argsort(all_spectra.fibermap["TARGETID"])
    nbpages = int(np.ceil((nspec_tile/nspecperfile)))
    for i_page in range(1,1+nbpages) :

        log.info(" * Page "+str(i_page)+" / "+str(nbpages))
        the_indices = sort_indices[(i_page-1)*nspecperfile:i_page*nspecperfile]            
        thespec = myspecselect.myspecselect(all_spectra, indices=the_indices)
        titlepage = titlepage_prefix+"_"+str(i_page)
        plotframes.plotspectra(thespec, with_noise=True, with_coaddcam=True, is_coadded=True, zcatalog=zcat,
                    title=titlepage, html_dir=html_dir, mask_type='CMX_TARGET', with_thumb_only_page=True)
    nspec_done += nspec_tile
        
    return nspec_done



def main(args) :
    
    log = get_logger()
    webdir = args.webdir
    assert args.frametype in ['sframe', 'cframe']
    if ( [args.exposure_list, args.exposure, args.tile, args.tile_list] ).count(None) != 3 :
        log.info("Specview_cmx_frames : Wrong set of input tiles/exposures. Exiting")
        return 0
    
    page_sorting = "tile"
    if (args.exposure_list is not None) or (args.exposure is not None) : 
        page_sorting = "exposure"
    else :
        if args.sort_exposures : page_sorting = "tile-expo"
        
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
        
    if page_sorting in ["tile-expo","tile"] :
        tile_subset = None
        if args.tile_list is not None :
            tile_subset = np.loadtxt(args.tile_list, dtype=str, comments='#')
        if args.tile is not None :
            tile_subset = [ args.tile ]
        merge_exposures = True if page_sorting=='tile' else False
        subset_db = tile_db(args.specprod_dir, frametype=args.frametype, tile_subset=tile_subset, merge_exposures=merge_exposures, petals=args.petals)
        if tile_subset is not None :
            tmplist = [ x['tile'] for x in subset_db ]
            missing_tiles = [ x for x in tile_subset if x not in tmplist ]
            for x in missing_tiles : log.info("Missing tile, cannot be processed : "+x)
        log.info(str(len(subset_db))+" tiles [exposures] to be processed")
            
            
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
        if page_sorting == 'tile-expo' :
            log.info("Working on tile "+the_subset['tile']+" (exposure "+the_subset['exposure']+")")
            fdir = os.path.join( args.specprod_dir, 'tiles', the_subset['tile'], the_subset['night'] )
            html_dir = os.path.join(webdir,"tiles_all",the_subset['tile'])
            titlepage_prefix = "tile"+the_subset['tile']+"_expo"+the_subset['exposure']
            if args.mask != None :
                html_dir = os.path.join(webdir, "tiles_"+args.mask, the_subset['tile'])
                titlepage_prefix = args.mask+"_"+titlepage_prefix
        if page_sorting == 'tile' :
            log.info("Working on tile "+the_subset['tile'])
            fdir = os.path.join( args.specprod_dir, 'tiles', the_subset['tile'] )
            html_dir = os.path.join(webdir,"tiles_all",the_subset['tile'])
            titlepage_prefix = "tile"+the_subset['tile']
            if args.mask != None :
                html_dir = os.path.join(webdir, "tiles_"+args.mask, the_subset['tile'])
                titlepage_prefix = args.mask+"_"+titlepage_prefix

        if not os.path.exists(html_dir) : 
            os.makedirs(html_dir)
        
        if page_sorting == 'tile' :
            nspec_added = page_subset_tile(fdir, the_subset, args.frametype, html_dir, titlepage_prefix, args.mask, log, args.nspecperfile, args.snrcut, with_zcatalog=args.with_zcatalog)
        else :
            nspec_added = page_subset_expo(fdir, the_subset['exposure'], args.frametype, the_subset['petals'], html_dir, titlepage_prefix, args.mask, log, args.nspecperfile, args.snrcut)
                    
        # Stop running if needed, only once a full exposure is completed
        nspec_done += nspec_added
        if args.nmax_spectra is not None :
            if nspec_done >= args.nmax_spectra :
                log.info(str(nspec_done)+" spectra done : no other exposure will be processed")
                break

        return 0



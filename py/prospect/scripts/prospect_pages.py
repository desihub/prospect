# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
====================================
prospect.scripts.specview_cmx_coadds
====================================

Create static html files to visually inspect DESI spectra.
"""

import os
import argparse
import numpy as np
from astropy.table import Table, vstack
import astropy.io.fits

import desispec.io
import desispec.spectra
from desiutil.log import get_logger

from ..viewer import plotspectra
from ..utilities import (create_targetdb, create_subsetdb, match_catalog_to_spectra, match_rrdetails_to_spectra,
                        load_spectra_zcat_from_targets, metadata_selection, get_subset_label, frames2spectra)

def _parse():

    parser = argparse.ArgumentParser(description='Create static html files to visually inspect DESI spectra.',
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    #- Single file input
    parser.add_argument('--spectra_files', help='[Mode: Explicit input files] Absolute path of file(s) with DESI spectra. All input spectra files must have exactly the same format (fibermap, extra, scores...). Frames are not supported.', nargs='+', type=str, default=None)
    parser.add_argument('--zcat_files', help='[Mode: Explicit input files] Absolute path of redshift catalog fits file(s), matched one-by-one to spectra_files', nargs='+', type=str, default=None)
    parser.add_argument('--redrock_details_files', help='[Mode: Explicit input files] Absolute path of detailed redrock file(s) (.h5), matched one-by-one to spectra_files', nargs='+', type=str, default=None)
    
    #- "Multi" file input (can select some data subsets based on tiles, expids...)
    parser.add_argument('--datadir', help='[Mode: Scan directory tree] Location of input directory tree', type=str, default=None)
    parser.add_argument('--dirtree_type', help='''[Mode: Scan directory tree] The following directory tree categories are supported:
                dirtree_type='pernight': {datadir}/{tileid}/{night}/{spectra_type}-{petal}-{tile}-{night}.fits ;
                dirtree_type='perexp': {datadir}/{tileid}/{expid}/{spectra_type}-{petal}-{tile}-exp{expid}.fits ;
                dirtree_type='cumulative': {datadir}/{tileid}/{night}/{spectra_type}-{petal}-{tile}-thru{night}.fits ;
                dirtree_type='exposures': {datadir}/{night}/{expid}/{spectra_type}-{band}{petal}-{expid}.fits ;
            Note that 'perexp' and 'exposures' are different.
            To use blanc/cascades 'all' (resp 'deep') coadds, use dirtree_type='pernight' and nights=['all'] (resp ['deep'])''', type=str, default=None)
    parser.add_argument('--spectra_type', help='[Mode: Scan directory tree] Prefix of spectra files: frame, cframe, sframe, coadd, spectra. [c/s]frames are needed when dirtree_type=exposures.', type=str, default='coadd')
    parser_tiles = parser.add_mutually_exclusive_group()
    parser_tiles.add_argument('--tiles', help='[Mode: Scan directory tree] Filter over tile[s]', nargs='+', type=str, default=None)
    parser_tiles.add_argument('--tile_list', help='[Mode: Scan directory tree] Filter over tile[s]: ASCII file with list of tiles, one per row', type=str, default=None)
    parser_nights= parser.add_mutually_exclusive_group()
    parser_nights.add_argument('--nights', help='[Mode: Scan directory tree] Filter over night[s]', nargs='+', type=str, default=None)
    parser_nights.add_argument('--night_list', help='[Mode: Scan directory tree] Filter over night[s]: ASCII file with list of nights, one per row', type=str, default=None)
    parser_expids = parser.add_mutually_exclusive_group()
    parser_expids.add_argument('--expids', help='[Mode: Scan directory tree] Filter over exposure[s]', nargs='+', type=str, default=None)
    parser_expids.add_argument('--expid_list', help='[Mode: Scan directory tree] Filter over exposure[s]: ASCII file with list of exposure, one per row', type=str, default=None)
    parser.add_argument('--petals', help='[Mode: Scan directory tree] Filter over petal[s], labelled 0 to 9', nargs='+', type=str, default=None)
    parser.add_argument('--with_zcatalog', help='[Mode: Scan directory tree] Include redshift fit result (from "redrock"/"zbest" fits files)', action='store_true')
    parser.add_argument('--with_multiple_models', help='[Mode: Scan directory tree] Display several models, using detailed redrock files (.h5)', action='store_true')
    parser.add_argument('--nmax_spectra', help='[Mode: Scan directory tree] Stop creating html pages once nmax_spectra are done', type=int, default=None)

    #- Generic parameters. The only mandatory one is outputdir
    parser.add_argument('--outputdir', '-o', help='Directory where output html pages are writen', type=str)
    parser.add_argument('--titlepage_prefix', help='Prefix for html page title', type=str, default='DESI_spectra')
    parser.add_argument('--nspecperfile', help='Number of spectra in each html page', type=int, default=50)
    parser.add_argument('--top_metadata', help="""List of fibermap's metadata to be highlighted (display in the most visible table).
        Note: if fibermap contains FIRST/LAST/NUM_XX, then including XX in top_metadata will display all of FIRST/LAST/NUM_XX.""", nargs='+', type=str, default=None)
    parser.add_argument('--vi_countdown', help='Countdown widget (in minutes)', type=int, default=-1)
    parser.add_argument('--with_thumbnail_only_pages', help='Create independent html pages with a thumbnail gallery, in addition to the main prospect pages. These additional pages are much faster to download than the main pages.', action='store_true')
    parser.add_argument('--mask_type', help='Targeting mask type', type=str, default='DESI_TARGET')
    parser.add_argument('--template_dir', help='Redrock template directory', type=str, default=None)
    parser.add_argument('--clean_fiberstatus', dest='clean_fiberstatus', help='Filter out spectra with FIBERSTATUS!=0 (even if a target list is provided)', action='store_true')
    parser.add_argument('--no-clean_fiberstatus', dest='clean_fiberstatus', action='store_false')
    parser.set_defaults(clean_fiberstatus=True)
    parser.add_argument('--clean_bad_fibers_cmx', help='[Specific to CMX conditions] Remove list of known bad fibers at CMX time.', action='store_true')

    #- Filtering at the spectra level
    parser.add_argument('--targeting_mask', help='Filter objects with a given targeting mask.', type=str, default=None)
    parser.add_argument('--snr_min', help='Filter objects with scores.MEDIAN_COADD_SNR_B+R+Z (or MEDIAN_CALIB_SNR_B+R+Z) in range [min, max]', type=float, default=None)
    parser.add_argument('--snr_max', type=float, default=None)
    parser.add_argument('--gmag_min', help='Filter objects with [dereddened] g-mag in range [min, max]', type=float, default=None)
    parser.add_argument('--gmag_max', type=float, default=None)
    parser.add_argument('--rmag_min', help='Filter objects with [dereddened] r-mag in range [min, max]', type=float, default=None)
    parser.add_argument('--rmag_max', type=float, default=None)
    parser.add_argument('--chi2_min', help='Filter objects with Delta_chi2 (from redshift catalog file) in range [min, max]', type=float, default=None)
    parser.add_argument('--chi2_max', type=float, default=None)
    parse_target = parser.add_mutually_exclusive_group()
    parse_target.add_argument('--targets', help='Filter over TARGETID', nargs='+', type=int, default=None)
    parse_target.add_argument('--target_list', help='''Filter over TARGETID: ASCII or FITS file providing the list of TARGETIDs. File format:
        - ASCII: simple list of TARGETIDs, one per row.
        - FITS: file whose first HDU is a table including a column named TARGETID.''', type=str, default=None)

    args = parser.parse_args()
    return args


def _filter_list(args, filter_name):
    """ Make list of 'filter_name'
        from either args.filter_names (values are explicitely given to parser) 
        or args.filter_name_list (read file)
    """
    list_output = None
    if vars(args)[filter_name+'s'] is not None:
        list_output = vars(args)[filter_name+'s']
    elif vars(args)[filter_name+'_list'] is not None:
        fname = vars(args)[filter_name+'_list']
        #- TARGETIDs can be read from a FITS file
        if (filter_name == 'target') and (fname[-5:] in ['.fits', '.FITS']):
            hlist = astropy.io.fits.open(fname)
            list_output = hlist[1].data['TARGETID']
        else:
            list_output = np.loadtxt(fname, dtype=str, comments='#')

    #- TARGETIDs should really be int64
    if filter_name == 'target' and list_output is not None:
        list_output = np.array(list_output, dtype='int64')

    return list_output


# List of bad fibers in CMX data (see eg SB / KD emails 23-24/03/2020)
_bad_fibers_cmx = [
    [1000, 1250], # r2 amp A glowing spot CCD defect
    [2250, 2500], # r4 amp D serial register CTE problem
    [4500, 4750]  # b9 amp C readout salt-and-pepper noise
]


def load_spectra_zcat_from_dbentry(db_entry, args, log, with_redrock_version=True):
    """ Load (spectra, zcat, redrock_cat) from a subset of DESI data (a subdirectory).
        Somewhat similar to `..utilities.load_spectra_zcat_from_targets()`, but this function
        is intended to be used only by `prospect_pages` script.
    """
    spectra, zcat, redrock_cat = None, None, None
    ztables, rrtables = [], []
    subset_label = get_subset_label(db_entry['subset'], args.dirtree_type)

    for petal in db_entry['petals']:
        the_dir = os.path.join(args.datadir, db_entry['dataset'], db_entry['subset'])
        if args.dirtree_type == 'exposures':
            #- Read frames and convert to Spectra
            file_labels = [ band+petal+'-'+db_entry['subset'] for band in ['b','r','z'] ]
            fnames = [ os.path.join(the_dir, args.spectra_type+"-"+label+".fits") for label in file_labels ]
            frames = [ desispec.io.read_frame(fname) for fname in fnames ]
            the_spec = frames2spectra(frames, with_scores=True)
        else:
            #- Read a single Spectra file directly
            file_label = '-'.join([petal, db_entry['dataset'], subset_label])
            the_spec = desispec.io.read_spectra(os.path.join(the_dir, args.spectra_type+"-"+file_label+".fits"))
        if args.with_zcatalog:
            if os.path.isfile(os.path.join(the_dir, "redrock-"+file_label+".fits")):
                redrock_is_pre_everest = False
                the_zcat = Table.read(os.path.join(the_dir, "redrock-"+file_label+".fits"), 'REDSHIFTS')
            else: # pre-everest Redrock file nomenclature
                redrock_is_pre_everest = True
                the_zcat = Table.read(os.path.join(the_dir, "zbest-"+file_label+".fits"), 'ZBEST')
            if with_redrock_version:
                if redrock_is_pre_everest:
                    hdulist = astropy.io.fits.open(os.path.join(the_dir, "zbest-"+file_label+".fits"))
                else:
                    hdulist = astropy.io.fits.open(os.path.join(the_dir, "redrock-"+file_label+".fits"))
                the_zcat['RRVER'] = hdulist[hdulist.index_of('PRIMARY')].header['RRVER']
        else:
            the_zcat = None
        #- Filtering (done after zcat is loaded, in order to include chi2_min/max !)
        the_spec = metadata_selection(the_spec, log=log, mask=args.targeting_mask, mask_type=args.mask_type,
                                      snr_range=[args.snr_min, args.snr_max], gmag_range=[args.gmag_min, args.gmag_max],
                                      rmag_range=[args.rmag_min, args.rmag_max], clean_fiberstatus=args.clean_fiberstatus,
                                      chi2_range=[args.chi2_min, args.chi2_max], zcat=the_zcat) # TODO dirty_mask_merge?
        #- Fiber-based filter in CMX data kept here for record only
        if args.clean_bad_fibers_cmx:
            fibers = np.arange(5000)
            for cut_fiber_range in _bad_fibers_cmx :
                fibers = fibers[ ( (fibers < cut_fiber_range[0]) | (fibers > cut_fiber_range[1]) )]
            try:
                the_spec = the_spec.select(fibers=fibers)
            except RuntimeError as select_err:
                the_spec = None
        if the_spec is None:
            continue
        #- Match zcat after spectra are filtered !
        if args.with_zcatalog:
            the_zcat = match_catalog_to_spectra(the_zcat, the_spec)
            ztables.append(the_zcat)
        if args.with_multiple_models:
            if redrock_is_pre_everest:
                rrfile = os.path.join(the_dir, "redrock-"+file_label+".h5")
            else:
                rrfile = os.path.join(the_dir, "rrdetails-"+file_label+".h5")
            the_rrcat = match_rrdetails_to_spectra(rrfile, the_spec, Nfit=None) # TODO Nfit?
            rrtables.append(the_rrcat)
        if spectra is None:
            spectra = the_spec
        elif the_spec is not None:
            #- Use update() instead of stack(): handle the case when different files'fibermaps differ.
            spectra.update(the_spec)

    if args.with_zcatalog:
        zcat = vstack(ztables)
    if args.with_multiple_models:
        redrock_cat = vstack(rrtables)
        
    return (spectra, zcat, redrock_cat)


def page_subset(spectra, nspecperfile, titlepage_prefix, viewer_params, log,
                zcat=None, redrock_cat=None, sort_by_targetid=False):
    """ Make a set of prospect html pages from (spectra, zcat, redrock_cat)
    """
    nspec_tot = spectra.num_spectra()
    log.info(" * Total nb spectra for this set of VI pages: "+str(nspec_tot))
    if sort_by_targetid:
        sort_indices = np.argsort(spectra.fibermap["TARGETID"])
    else:
        sort_indices = np.arange(nspec_tot)
    nbpages = int(np.ceil((nspec_tot/nspecperfile)))
    for i_page in range(1,1+nbpages) :
        log.info(" * Page "+str(i_page)+" / "+str(nbpages))
        the_indices = sort_indices[(i_page-1)*nspecperfile:i_page*nspecperfile]
        thespec = spectra[the_indices]
        if zcat is not None:
            the_zcat = match_catalog_to_spectra(zcat, thespec)
        else:
            the_zcat = None
        if redrock_cat is not None :
            the_rrtable = match_catalog_to_spectra(redrock_cat, thespec)
        else :
            the_rrtable = None
        titlepage = titlepage_prefix
        if nbpages>1: titlepage += ("_"+str(i_page))
        plotspectra(thespec, zcatalog=the_zcat, redrock_cat=the_rrtable,
                    title=titlepage, html_dir=viewer_params['html_dir'],
                    mask_type=viewer_params['mask_type'], top_metadata=viewer_params['top_metadata'],
                    template_dir=viewer_params['template_dir'], num_approx_fits=viewer_params['num_approx_fits'],
                    with_full_2ndfit=viewer_params['with_full_2ndfit'], vi_countdown=viewer_params['vi_countdown'],
                    with_thumb_only_page=viewer_params['with_thumb_only_page'])

    return nspec_tot


def main():
    args = _parse()
    log = get_logger()
    
    #- Two ways to provide input files
    if args.spectra_files is None:
        input_mode = 'scan-dirtree'
    else:
        input_mode = 'explicit-files'

    if (args.targets is None) and (args.target_list is None):
        target_filter = False
    else:
        target_filter = True
        #- If a list of TARGETIDs is given, one should not ask any other filtering at the spectra level
        for the_arg in ['targeting_mask', 'snr_min', 'snr_max', 'gmag_min', 'gmag_max', 'rmag_min', 'rmag_max', 'chi2_min', 'chi2_max', 'nmax_spectra']:
            if vars(args)[the_arg] is not None:
                raise ValueError('Argument(s) not allowed with --targets/--target_list: '+the_arg)

    #- Parameters provided as is to viewer.plotspectra:
    viewer_params = {
        'html_dir': args.outputdir,
        'mask_type': args.mask_type,
        'top_metadata': args.top_metadata,
        'template_dir': args.template_dir,
        'vi_countdown': args.vi_countdown,
        'num_approx_fits': None,
        'with_full_2ndfit': False,
        'with_thumb_only_page': args.with_thumbnail_only_pages,
    }


    #############################
    #- Mode: Explicit input files
    #############################
    if input_mode == 'explicit-files':
        for the_arg in ['datadir', 'dirtree_type', 'tiles', 'tile_list',
                        'nights', 'night_list', 'expids', 'expid_list', 'petals', 'nmax_spectra']:
            if vars(args)[the_arg] is not None:
                raise ValueError('Argument not allowed in mode "Explicit input files": '+the_arg)
        if (args.redrock_details_files is not None) and (args.zcat_files is None):
            raise ValueError('Argument `zcat_files` is needed if `redrock_details_files` is set')
        n_specfiles = len(args.spectra_files)
        if args.zcat_files is not None :
            if len(args.zcat_files)!=n_specfiles:
                raise ValueError('Number of zcat_files does not match number of input spectra_files')
        if args.redrock_details_files is not None :
            viewer_params['num_approx_fits'] = 4 # TODO un-hardcode ?
            viewer_params['with_full_2ndfit'] = True # TODO un-hardcode ?
            if len(args.redrock_details_files)!=n_specfiles:
                raise ValueError('Number of redrock_details_files does not match number of input spectra_files')
 
        log.info('Prospect_pages: start reading data [mode: Explicit input files]')
        #- Read input file(s)
        spectra_list, zcat_list, redrock_list = [], [], []
        for i_file in range(n_specfiles):
            spectra = desispec.io.read_spectra(args.spectra_files[i_file])
            spectra_list.append(spectra)
            if args.zcat_files is not None:
                try:
                    zcat = Table.read(args.zcat_files[i_file],'REDSHIFTS')
                except KeyError: # pre-everest Redrock file nomenclature
                    zcat = Table.read(args.zcat_files[i_file],'ZBEST')
                #- Add redrock version to zcat
                hdulist = astropy.io.fits.open(args.zcat_files[i_file])
                zcat['RRVER'] = hdulist[hdulist.index_of('PRIMARY')].header['RRVER']
                zcat_list.append(zcat)
            if args.redrock_details_files is not None:
                redrock_cat = match_rrdetails_to_spectra(args.redrock_details_files[i_file], spectra)
                redrock_list.append(redrock_cat)
        spectra = desispec.spectra.stack(spectra_list)
        if args.zcat_files is not None:
            zcat = vstack(zcat_list)
        else:
            zcat = None
        if args.redrock_details_files is not None:
            redrock_cat = vstack(redrock_list)
        else:
            redrock_cat = None
        #- Filtering: generic metadata
        spectra = metadata_selection(spectra, log=log, mask=args.targeting_mask, mask_type=args.mask_type,
                                     snr_range=[args.snr_min, args.snr_max], gmag_range=[args.gmag_min, args.gmag_max],
                                     rmag_range=[args.rmag_min, args.rmag_max], clean_fiberstatus=args.clean_fiberstatus,
                                     chi2_range=[args.chi2_min, args.chi2_max], zcat=zcat) # TODO dirty_mask_merge?
        if spectra is None:
            log.info("Prospect_pages: no spectra after filtering criteria -> End of task.")
            return 0
        #- Filtering: list of TARGETIDs. NB Here spectra are just filtered, ie. not reordered according to input targets
        if target_filter:
            targetids = _filter_list(args, 'target')
            spectra = spectra.select(targets=targetids)
        #- Run viewer.plotspectra()
        log.info('Prospect_pages: start creating html page(s)')
        n_done = page_subset(spectra, args.nspecperfile, args.titlepage_prefix, viewer_params, log, zcat=zcat, redrock_cat=redrock_cat)


    ############################
    #- Mode: Scan directory tree
    ############################
    if input_mode == 'scan-dirtree':
        if any([ x is not None for x in [args.spectra_files, args.zcat_files, args.redrock_details_files] ]):
            raise ValueError('Argument not allowed in "Scan directory tree" mode: spectra/zcat/redrock_details_files')
        if any([ x is None for x in [args.datadir, args.dirtree_type] ]):
            raise ValueError('Missing parameter in "Scan directory tree" mode: datadir/dirtree_type')
        if args.with_multiple_models and not args.with_zcatalog:
            raise ValueError('Need `with_zcatalog` if `with_multiple_models` is set')
        if target_filter and not args.with_zcatalog:
            raise ValueError('Need `with_zcatalog` if TARGETID-based filtering in multi-file input mode')
        if args.spectra_type in ['frame', 'cframe', 'sframe']:
            if (args.dirtree_type!='exposures') or args.with_zcatalog:
                raise ValueError('Argument(s) dirtree_type/with_zcatalog not consistent with frame-like spectra_type.')
        elif args.spectra_type in ['coadd', 'spectra']:
            if args.dirtree_type == 'exposures':
                raise ValueError('Argument dirtree_type not consistent with spectra_type.')
        else:
            raise ValueError('spectra_type not supported: '+args.spectra_type)

        tiles = _filter_list(args, 'tile')
        nights = _filter_list(args, 'night')
        expids = _filter_list(args, 'expid')
        subset_db = create_subsetdb(args.datadir, dirtree_type=args.dirtree_type, spectra_type=args.spectra_type,
                                   tiles=tiles, nights=nights, expids=expids, petals=args.petals,
                                   with_zcat=args.with_zcatalog)

        #- Spectra from a list of TARGETIDs
        if target_filter:
            log.info('Prospect_pages: start reading data [mode: Scan directory tree, with target selection]')
            target_db = create_targetdb(args.datadir, subset_db, dirtree_type=args.dirtree_type)
            targetids = _filter_list(args, 'target')
            if args.with_multiple_models:
                spectra, zcat, redrock_cat = load_spectra_zcat_from_targets(targetids, args.datadir, target_db,
                                                    dirtree_type=args.dirtree_type, with_redrock_details=True)
                viewer_params['num_approx_fits'] = 4 # TODO un-hardcode ?
                viewer_params['with_full_2ndfit'] = True # TODO un-hardcode ?
            else:
                spectra, zcat = load_spectra_zcat_from_targets(targetids, args.datadir, target_db,
                                                    dirtree_type=args.dirtree_type, with_redrock_details=False)
                redrock_cat = None
            log.info('Prospect_pages: start creating html page(s)')
            n_done = page_subset(spectra, args.nspecperfile, args.titlepage_prefix, viewer_params, log, zcat=zcat, redrock_cat=redrock_cat)
        
        #- All spectra, possibly filtered based on some metadata
        else:
            log.info('Prospect_pages: start reading data [mode: Scan directory tree]')
            n_done = 0
            for db_entry in subset_db:
                log.info("Tile "+db_entry['dataset']+" - subset "+db_entry['subset'])
                spectra, zcat, redrock_cat = load_spectra_zcat_from_dbentry(db_entry, args, log)
                if spectra is None:
                    log.info('No spectra found for this subset')
                    continue
                if redrock_cat is not None:
                    viewer_params['num_approx_fits'] = 4 # TODO un-hardcode ?
                    viewer_params['with_full_2ndfit'] = True # TODO un-hardcode ?
                #- Associate a subdirectory for each individual subset:
                html_subdir = db_entry['dataset']+'-'+get_subset_label(db_entry['subset'], args.dirtree_type)
                viewer_params['html_dir'] = os.path.join(args.outputdir, html_subdir)
                os.makedirs(viewer_params['html_dir'], exist_ok=True)
                titlepage_prefix = args.titlepage_prefix + '_' + html_subdir
                nspec = page_subset(spectra, args.nspecperfile, titlepage_prefix, viewer_params, log, zcat=zcat, redrock_cat=redrock_cat)
                n_done += nspec
                if args.nmax_spectra is not None:
                    if n_done >= args.nmax_spectra:
                        log.info('Prospect_pages: reached nmax_spectra -> End of task.')
                        return 0

    log.info('Prospect_pages: all done -> End of task.')
    return 0


    


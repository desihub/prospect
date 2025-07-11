# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
===============================
prospect.scripts.prospect_pages
===============================

Create static html files to visually inspect DESI spectra.
"""

import os
import argparse
import numpy as np
from astropy.table import Table, vstack
import astropy.io.fits

import desispec.io
import desispec.spectra
from desispec.io.util import healpix_subdirectory
from desiutil.log import get_logger

from ..viewer import plotspectra
from ..utilities import (create_targetdb, create_subsetdb, match_catalog_to_spectra, match_rrdetails_to_spectra, file_or_gz,
                        file_or_gz_exists, load_spectra_zcat_from_targets, metadata_selection, get_subset_label, frames2spectra)

def _parse():

    parser = argparse.ArgumentParser(description='Create static html files to visually inspect DESI spectra. Spectra are selected from '
                                     'either a list of files (Mode: Explicit input files), or a given DESI directory tree to be parsed by prospect '
                                     ' (Mode: Scan directory tree)', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    #- Explicit file input
    parser.add_argument('--spectra_files', help='[Mode: Explicit input files] Absolute path of file(s) with DESI spectra. All input spectra files must have exactly the same format (fibermap, extra, scores...). Frames are not supported.', nargs='+', type=str, default=None)
    parser.add_argument('--spectra_file_list', help='[Mode: Explicit input files] ASCII file with list of DESI spectra files, one per row', type=str, default=None)
    parser.add_argument('--zcat_files', help='[Mode: Explicit input files] Absolute path of redshift catalog fits file(s), matched one-by-one to spectra_files', nargs='+', type=str, default=None)
    parser.add_argument('--zcat_file_list', help='[Mode: Explicit input files] ASCII file with list of redshift catalog files, one per row', type=str, default=None)
    parser.add_argument('--redrock_details_files', help='[Mode: Explicit input files] Absolute path of detailed redrock file(s) (.h5), matched one-by-one to spectra_files', nargs='+', type=str, default=None)
    parser.add_argument('--redrock_details_file_list', help='[Mode: Explicit input files] ASCII file with list of detailed redrock files, one per row', type=str, default=None)

    #- Selection-based file input (select data subsets based on tiles, expids...)
    parser.add_argument('--datadir', help='[Mode: Scan directory tree] Location of input directory tree (eg. $DESI_SPECTRO_REDUX/iron/healpix)', type=str, default=None)
    parser.add_argument('--dirtree_type', help='''[Mode: Scan directory tree] The following directory tree categories are supported:
                dirtree_type='healpix': {datadir}/{survey}/{program}/{pixel//100}/{pixel}/{spectra_type}-{survey}-{program}-{pixel}.fits ;
                dirtree_type='pernight': {datadir}/{tileid}/{night}/{spectra_type}-{petal}-{tile}-{night}.fits ;
                dirtree_type='perexp': {datadir}/{tileid}/{expid}/{spectra_type}-{petal}-{tile}-exp{expid}.fits ;
                dirtree_type='cumulative': {datadir}/{tileid}/{night}/{spectra_type}-{petal}-{tile}-thru{night}.fits ;
                dirtree_type='exposures': {datadir}/{night}/{expid}/{spectra_type}-{band}{petal}-{expid}.fits ;
            Note that 'perexp' and 'exposures' are different. Fits and fits.gz extensions are ok.
            To use blanc/cascades 'all' (resp 'deep') coadds, use dirtree_type='pernight' and nights=['all'] (resp ['deep'])''', type=str, default=None)
    parser.add_argument('--spectra_type', help='[Mode: Scan directory tree] Prefix of spectra files: frame, cframe, sframe, coadd, spectra. [c/s]frames are needed when dirtree_type=exposures.', type=str, default='coadd')
    parser.add_argument('--survey_program', help='[Mode: Scan directory tree] Set survey and program (eg. --survey_program sv3 dark) when dirtree_type=healpix', nargs='+', type=str, default=None)
    parser_pixels = parser.add_mutually_exclusive_group()
    parser_pixels.add_argument('--pixels', help='[Mode: Scan directory tree] Filter over pixel[s]', nargs='+', type=str, default=None)
    parser_pixels.add_argument('--pixel_list', help='[Mode: Scan directory tree] Filter over pixel[s]: ASCII file with list of pixels, one per row', type=str, default=None)
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
    parser.add_argument('--nspec_per_page', help='Number of spectra in each html page', type=int, default=50)
    parser.add_argument('--with_thumbnail_only_pages', help='Create independent html pages with a thumbnail gallery, in addition to the main prospect pages. These additional pages are much faster to download than the main pages.', action='store_true')
    parser.add_argument('--mask_type', help='Targeting mask type', type=str, default='DESI_TARGET')
    parser.add_argument('--template_dir', help='Redrock template directory', type=str, default=None)
    parser.add_argument('--std_template_file', help='''File containing "standard templates" to display in viewer.
                Format: N templates in a single fits file, containing N tables (HDU 1 to N) exclusively.
                The table associated to template TMPLT should contain two columns, named 'wave_TMPLT' and 'flux_TMPLT'.
                Wavelength arrays should be regularly, linearly binned.
            Prospect default is in data/std_template.fits.''', type=str, default=None)
    parser.add_argument('--no_clean_fiberstatus', dest='clean_fiberstatus', help='Do not filter out spectra with FIBERSTATUS!=0', action='store_false')
    parser.add_argument('--select_bad_fiberstatus', help='[For debugging] Select spectra with FIBERSTATUS!=0', action='store_true')
    parser.add_argument('--clean_bad_fibers_cmx', help='[Specific to CMX conditions] Remove list of known bad fibers at CMX time.', action='store_true')

    #- VI page display parameters
    parser.add_argument('--titlepage_prefix', help='Prefix for html page title', type=str, default='DESI_spectra')
    parser.add_argument('--top_metadata', help="""List of fibermap's metadata to be highlighted (display in the most visible table).
        Note: if fibermap contains FIRST/LAST/NUM_XX, then including XX in top_metadata will display all of FIRST/LAST/NUM_XX.""", nargs='+', type=str, default=None)
    parser.add_argument('--colors', help="""Customize the curve's colors: 3 colors should be given, associated respectively to the coadded data, the model and the noise.""", nargs='+', type=str, default=None)
    parser.add_argument('--no_imaging', dest='with_imaging', help='Do not include thumb images from https://www.legacysurvey.org/viewer', action='store_false')
    parser.add_argument('--no_noise', dest='with_noise', help='Do not display noise vectors associated to spectra', action='store_false')
    parser.add_argument('--no_other_model', dest='with_other_model', help="""Do not display the 'other model' curve""", action='store_false')
    parser.add_argument('--no_thumb_tab', dest='with_thumb_tab', help='Do not include a tab with spectra thumbnails', action='store_false')
    parser.add_argument('--no_vi_widgets', dest='with_vi_widgets', help='Do not include widgets used to enter VI information', action='store_false')
    parser.add_argument('--no_coaddcam', dest='with_coaddcam', help='Do not include camera-coaddition (DESI only)', action='store_false')
    parser.add_argument('--vi_countdown', help='Countdown widget (in minutes)', type=int, default=-1)
    parser.add_argument('--no_full_2ndfit', dest='with_full_2ndfit', help='Compute and display the second best-fit model without approximation (when available)', action='store_false')
    parser.add_argument('--num_approx_fits', help='Number of best-fit models to display', type=int, default=4)
    parser.add_argument('--zmax_slider', help='Maximum range of the redshift slider widget', type=float, default=5.0)

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
        if args.dirtree_type == 'healpix':
            nside = 64 # dummy, currently
            the_dir = os.path.join(args.datadir, db_entry['dataset'][0], db_entry['dataset'][1],
                                   healpix_subdirectory(nside, int(db_entry['subset'])))
            file_label = '-'.join([db_entry['dataset'][0], db_entry['dataset'][1], subset_label])
        else:
            the_dir = os.path.join(args.datadir, db_entry['dataset'], db_entry['subset'])
            if args.dirtree_type == 'exposures':
                file_label = [ band+petal+'-'+db_entry['subset'] for band in ['b','r','z'] ]
            else:
                file_label = '-'.join([petal, db_entry['dataset'], subset_label])
        if args.dirtree_type == 'exposures':
            #- Read frames and convert to Spectra
            fnames = [ os.path.join(the_dir, args.spectra_type+"-"+label+".fits") for label in file_label ]
            frames = [ desispec.io.read_frame(fname) for fname in fnames ]
            the_spec = frames2spectra(frames, with_scores=True)
        else:
            #- Read a single Spectra file directly
            the_spec = desispec.io.read_spectra(os.path.join(the_dir, args.spectra_type+"-"+file_label+".fits"), single=True)
        if args.with_zcatalog:
            if file_or_gz_exists(os.path.join(the_dir, "redrock-"+file_label+".fits")):
                redrock_is_pre_everest = False
                the_zcat = Table.read(file_or_gz(os.path.join(the_dir, "redrock-"+file_label+".fits")), 'REDSHIFTS')
            else:  # pre-everest Redrock file nomenclature
                redrock_is_pre_everest = True
                the_zcat = Table.read(file_or_gz(os.path.join(the_dir, "zbest-"+file_label+".fits")), 'ZBEST')
            if hasattr(the_zcat['SUBTYPE'], 'mask'):  # work around Table auto-masking in astropy 5
                blanksubtype = the_zcat['SUBTYPE'].mask
                the_zcat['SUBTYPE'][blanksubtype] = ''
            if with_redrock_version:
                if redrock_is_pre_everest:
                    hdulist = astropy.io.fits.open(file_or_gz(os.path.join(the_dir, "zbest-"+file_label+".fits")))
                else:
                    hdulist = astropy.io.fits.open(file_or_gz(os.path.join(the_dir, "redrock-"+file_label+".fits")))
                the_zcat['RRVER'] = hdulist[hdulist.index_of('PRIMARY')].header['RRVER']
        else:
            the_zcat = None
        #- Filtering (done after zcat is loaded, in order to include chi2_min/max !)
        the_spec = metadata_selection(the_spec, log=log, mask=args.targeting_mask, mask_type=args.mask_type,
                            snr_range=[args.snr_min, args.snr_max], gmag_range=[args.gmag_min, args.gmag_max],
                            rmag_range=[args.rmag_min, args.rmag_max], chi2_range=[args.chi2_min, args.chi2_max],
                            clean_fiberstatus=args.clean_fiberstatus, select_bad_fiberstatus=args.select_bad_fiberstatus,
                            zcat=the_zcat) # TODO dirty_mask_merge?
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


def page_subset(spectra, nspec_per_page, titlepage_prefix, viewer_params, log,
                zcat=None, redrock_cat=None, sort_by_targetid=False):
    """ Make a set of prospect html pages from (spectra, zcat, redrock_cat)
        zcat and redrock_cat must be entry-matched to spectra
    """
    nspec_tot = spectra.num_spectra()
    log.info(" * Total nb spectra for this set of VI pages: "+str(nspec_tot))
    if sort_by_targetid:
        sort_indices = np.argsort(spectra.fibermap["TARGETID"])
    else:
        sort_indices = np.arange(nspec_tot)
    nbpages = int(np.ceil((nspec_tot/nspec_per_page)))
    for i_page in range(1,1+nbpages) :
        log.info(" * Page "+str(i_page)+" / "+str(nbpages))
        the_indices = sort_indices[(i_page-1)*nspec_per_page:i_page*nspec_per_page]
        thespec = spectra[the_indices]
        if zcat is not None:
            the_zcat = zcat[the_indices]
        else:
            the_zcat = None
        if redrock_cat is not None :
            the_rrtable = redrock_cat[the_indices]
        else :
            the_rrtable = None
        titlepage = titlepage_prefix
        if nbpages>1: titlepage += ("_"+str(i_page))
        plotspectra(thespec, zcatalog=the_zcat, redrock_cat=the_rrtable,
                    title=titlepage, **viewer_params)

    return nspec_tot


def main():
    """Entry-point for :command:`prospect_pages`.
    """
    args = _parse()
    log = get_logger()

    #- Two ways to provide input files
    if (args.spectra_files is None) and (args.spectra_file_list is None):
        input_mode = 'scan-dirtree'
    else:
        input_mode = 'explicit-files'

    if (args.targets is None) and (args.target_list is None):
        target_filter = False
    else:
        target_filter = True
        targetids = _filter_list(args, 'target')
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
        'num_approx_fits': args.num_approx_fits,
        'with_full_2ndfit': args.with_full_2ndfit,
        'with_thumb_only_page': args.with_thumbnail_only_pages,
        'std_template_file': args.std_template_file,
        'colors': args.colors,
        'with_imaging': args.with_imaging,
        'with_noise': args.with_noise,
        'with_other_model': args.with_other_model,
        'with_thumb_tab': args.with_thumb_tab,
        'with_vi_widgets': args.with_vi_widgets,
        'with_coaddcam': args.with_coaddcam,
        'zmax_slider': args.zmax_slider
    }


    #############################
    #- Mode: Explicit input files
    #############################
    if input_mode == 'explicit-files':
        for the_arg in ['datadir', 'dirtree_type', 'tiles', 'tile_list',
                        'nights', 'night_list', 'expids', 'expid_list', 'petals', 'nmax_spectra']:
            if vars(args)[the_arg] is not None:
                raise ValueError('Argument not allowed in mode "Explicit input files": '+the_arg)
        spectra_files = _filter_list(args, 'spectra_file')
        zcat_files = _filter_list(args, 'zcat_file')
        redrock_details_files = _filter_list(args, 'redrock_details_file')
        if (redrock_details_files is not None) and (zcat_files is None):
            raise ValueError('Argument `zcat_files` is needed if `redrock_details_files` is set')
        n_specfiles = len(spectra_files)
        if zcat_files is not None :
            if len(zcat_files)!=n_specfiles:
                raise ValueError('Number of zcat_files does not match number of input spectra_files')
        if redrock_details_files is not None :
            if len(redrock_details_files)!=n_specfiles:
                raise ValueError('Number of redrock_details_files does not match number of input spectra_files')

        log.info('Prospect_pages: start reading data [mode: Explicit input files]')
        #- Read input file(s)
        spectra_list, zcat_list, redrock_list = [], [], []
        for i_file in range(n_specfiles):
            spectra = desispec.io.read_spectra(spectra_files[i_file], single=True)
            if zcat_files is not None:
                try:
                    zcat = Table.read(zcat_files[i_file], 'REDSHIFTS')
                except KeyError as e:  # pre-everest Redrock file nomenclature
                    zcat = Table.read(zcat_files[i_file], 'ZBEST')
                if hasattr(zcat['SUBTYPE'], 'mask'):  # work around Table auto-masking in astropy 5
                    blanksubtype = zcat['SUBTYPE'].mask
                    zcat['SUBTYPE'][blanksubtype] = ''
                #- Add redrock version to zcat
                hdulist = astropy.io.fits.open(zcat_files[i_file])
                zcat['RRVER'] = hdulist[hdulist.index_of('PRIMARY')].header['RRVER']
            else:
                zcat = None
            #- Filtering: generic metadata
            spectra, indx = metadata_selection(spectra, log=log, mask=args.targeting_mask, mask_type=args.mask_type,
                            snr_range=[args.snr_min, args.snr_max], gmag_range=[args.gmag_min, args.gmag_max],
                            rmag_range=[args.rmag_min, args.rmag_max], chi2_range=[args.chi2_min, args.chi2_max],
                            clean_fiberstatus=args.clean_fiberstatus, select_bad_fiberstatus=args.select_bad_fiberstatus,
                            zcat=zcat, return_index=True)
            if zcat_files is not None:
                zcat = zcat[indx]
            #- Filtering: targetids
            # NB Here spectra are just filtered, ie. not reordered according to input targets
            if target_filter:
                try:
                    spectra, indx = spectra.select(targets=targetids, return_index=True)
                except RuntimeError:   # this happens if no TARGETID is matched:
                    spectra, indx = None, np.array([], dtype='int64')
                if zcat_files is not None:
                    zcat = zcat[indx]
            if spectra is not None:
                spectra_list.append(spectra)
                if zcat_files is not None:
                    zcat_list.append(zcat)
                if redrock_details_files is not None:
                    redrock_cat = match_rrdetails_to_spectra(redrock_details_files[i_file], spectra)
                    redrock_list.append(redrock_cat)
        if len(spectra_list) == 0:
            log.info("Prospect_pages: no spectra after filtering criteria -> End of task.")
            return 0
        spectra = desispec.spectra.stack(spectra_list)
        del spectra_list[:]   # should roughly divide memory usage by two
        if zcat_files is not None:
            zcat = vstack(zcat_list)
        else:
            zcat = None
        if redrock_details_files is not None:
            redrock_cat = vstack(redrock_list)
        else:
            redrock_cat = None

        #- Run viewer.plotspectra()
        log.info('Prospect_pages: start creating html page(s)')
        n_done = page_subset(spectra, args.nspec_per_page, args.titlepage_prefix, viewer_params, log, zcat=zcat, redrock_cat=redrock_cat)


    ############################
    #- Mode: Scan directory tree
    ############################
    if input_mode == 'scan-dirtree':
        if any([ x is not None for x in [args.spectra_files, args.zcat_files, args.redrock_details_files,
                            args.spectra_file_list, args.zcat_file_list, args.redrock_details_file_list] ]):
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
        pixels = _filter_list(args, 'pixel')
        subset_db = create_subsetdb(args.datadir, dirtree_type=args.dirtree_type, spectra_type=args.spectra_type,
                                   tiles=tiles, nights=nights, expids=expids, pixels=pixels, petals=args.petals,
                                   survey_program=args.survey_program, with_zcat=args.with_zcatalog)

        #- Spectra from a list of TARGETIDs
        if target_filter:
            log.info('Prospect_pages: start reading data [mode: Scan directory tree, with target selection]')
            target_db = create_targetdb(args.datadir, subset_db, dirtree_type=args.dirtree_type)
            if args.with_multiple_models:
                spectra, zcat, redrock_cat = load_spectra_zcat_from_targets(targetids, args.datadir, target_db,
                                                    dirtree_type=args.dirtree_type, with_redrock_details=True)
            else:
                spectra, zcat = load_spectra_zcat_from_targets(targetids, args.datadir, target_db,
                                                    dirtree_type=args.dirtree_type, with_redrock_details=False)
                redrock_cat = None
            log.info('Prospect_pages: start creating html page(s)')
            n_done = page_subset(spectra, args.nspec_per_page, args.titlepage_prefix, viewer_params, log, zcat=zcat, redrock_cat=redrock_cat)

        #- All spectra, possibly filtered based on some metadata
        else:
            log.info('Prospect_pages: start reading data [mode: Scan directory tree]')
            n_done = 0
            for db_entry in subset_db:
                if args.dirtree_type == 'healpix':
                    dataset = '-'.join(db_entry['dataset'])
                else:
                    dataset = db_entry['dataset']
                log.info("Dataset "+dataset+" - subset "+db_entry['subset'])
                spectra, zcat, redrock_cat = load_spectra_zcat_from_dbentry(db_entry, args, log)
                if spectra is None:
                    log.info('No spectra found for this subset')
                    continue
                #- Associate a subdirectory for each individual subset:
                html_subdir = dataset+'-'+get_subset_label(db_entry['subset'], args.dirtree_type)
                viewer_params['html_dir'] = os.path.join(args.outputdir, html_subdir)
                os.makedirs(viewer_params['html_dir'], exist_ok=True)
                titlepage_prefix = args.titlepage_prefix + '_' + html_subdir
                nspec = page_subset(spectra, args.nspec_per_page, titlepage_prefix, viewer_params, log, zcat=zcat, redrock_cat=redrock_cat)
                n_done += nspec
                if args.nmax_spectra is not None:
                    if n_done >= args.nmax_spectra:
                        log.info('Prospect_pages: reached nmax_spectra -> End of task.')
                        return 0

    log.info('Prospect_pages: all done -> End of task.')
    return 0





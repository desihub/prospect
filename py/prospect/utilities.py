# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
==================
prospect.utilities
==================

Utility functions for prospect.
"""

import os, glob
from pkg_resources import resource_string, resource_listdir

import numpy as np
import astropy.io.fits
from astropy.table import Table, vstack
import scipy.ndimage.filters

_desiutil_imported = True
try:
    from desiutil.log import get_logger
except ImportError:
    _desiutil_imported = False

_desispec_imported = True
try:
    import desispec.spectra
    import desispec.frame
except ImportError:
    _desispec_imported = False

_desitarget_imported = True
try:
    from desitarget.targetmask import desi_mask, bgs_mask, mws_mask, scnd_mask
    from desitarget.cmx.cmx_targetmask import cmx_mask
    from desitarget.sv1.sv1_targetmask import desi_mask as sv1_desi_mask
    from desitarget.sv1.sv1_targetmask import bgs_mask as sv1_bgs_mask
    from desitarget.sv1.sv1_targetmask import mws_mask as sv1_mws_mask
    from desitarget.sv1.sv1_targetmask import scnd_mask as sv1_scnd_mask
    from desitarget.sv2.sv2_targetmask import desi_mask as sv2_desi_mask
    from desitarget.sv2.sv2_targetmask import bgs_mask as sv2_bgs_mask
    from desitarget.sv2.sv2_targetmask import mws_mask as sv2_mws_mask
    from desitarget.sv2.sv2_targetmask import scnd_mask as sv2_scnd_mask
    from desitarget.sv3.sv3_targetmask import desi_mask as sv3_desi_mask
    from desitarget.sv3.sv3_targetmask import bgs_mask as sv3_bgs_mask
    from desitarget.sv3.sv3_targetmask import mws_mask as sv3_mws_mask
    from desitarget.sv3.sv3_targetmask import scnd_mask as sv3_scnd_mask
    supported_desitarget_masks = {
        'DESI_TARGET': desi_mask,
        'BGS_TARGET': bgs_mask,
        'MWS_TARGET': mws_mask,
        'SECONDARY_TARGET': scnd_mask,
        'CMX_TARGET': cmx_mask,
        'SV1_DESI_TARGET': sv1_desi_mask,
        'SV1_BGS_TARGET': sv1_bgs_mask,
        'SV1_MWS_TARGET': sv1_mws_mask,
        'SV1_SCND_TARGET': sv1_scnd_mask,
        'SV2_DESI_TARGET': sv2_desi_mask,
        'SV2_BGS_TARGET': sv2_bgs_mask,
        'SV2_MWS_TARGET': sv2_mws_mask,
        'SV2_SCND_TARGET': sv2_scnd_mask,
        'SV3_DESI_TARGET': sv3_desi_mask,
        'SV3_BGS_TARGET': sv3_bgs_mask,
        'SV3_MWS_TARGET': sv3_mws_mask,
        'SV3_SCND_TARGET': sv3_scnd_mask,
        }
except ImportError:
    _desitarget_imported = False
    supported_desitarget_masks = dict()

_redrock_imported = True
try:
    import redrock.results
except ImportError:
    _redrock_imported = False

vi_flags = [
    # Definition of VI flags
    # shortlabels for "issue" flags must be a unique single-letter identifier
    {"label" : "4", "type" : "quality", "description" : "Confident classification: two or more secure features."},
    {"label" : "3", "type" : "quality", "description" : "Probable classification: at least one secure spectral feature + continuum or many weak spectral features."},
    {"label" : "2", "type" : "quality", "description" : "Possible classification: one strong spectral feature but unsure what it is."},
    {"label" : "1", "type" : "quality", "description" : "Unlikely classification: clear signal but features are unidentified."},
    {"label" : "0", "type" : "quality", "description" : "Nothing there, no signal."},
    {"label" : "Bad redshift fit", "shortlabel" : "R", "type" : "issue", "description" : "Mis-estimation of redshift by the pipeline fitter"},
    {"label" : "Bad spectype fit", "shortlabel" : "C", "type" : "issue", "description" : "Mis-identification of spectral type from the best-fit pipeline solution; e.g., star vs QSO..."},
    {"label" : "Bad spectrum", "shortlabel" : "S", "type" : "issue", "description" : "Bad spectrum; e.g. strong cosmic/skyline subtraction residuals."}
]

vi_file_fields = [
    # Contents of VI files: [
    #      field name (in VI file header),
    #      associated variable in viewer_cds.cds_metadata,
    #      dtype in VI file
    #      default value ]
    # Ordered list
    ["TARGETID", "TARGETID", "i8", -1],
    ["EXPID", "EXPID", "i4", -1],
    ["NIGHT", "NIGHT", "i4", -1],
    ["TILEID", "TILEID", "i4", -1],
    ["Spec_version", "spec_version", "U16", "-1"],
    ["Redrock_version", "redrock_version", "U16", "-1"],
    ["Template_version", "template_version", "U16", "-1"],
    ["Redrock_spectype", "SPECTYPE", "U10", ""],
    ["Redrock_z", "Z", "U6", "-1"],
    ["Redrock_deltachi2", "DELTACHI2", "U10", "-1"],
    ["VI_scanner", "VI_scanner", "U10", " "],
    ["VI_quality", "VI_quality_flag", "U2", "-1"],
    ["VI_issue", "VI_issue_flag", "U3", ""],
    ["VI_z", "VI_z", "U6", ""],
    ["VI_spectype", "VI_spectype", "U10", ""],
    ["VI_comment", "VI_comment", "U100", ""]
]

vi_spectypes =[
    # List of spectral types to fill in VI categories
    # in principle, it should match somehow redrock spectypes...
    "STAR",
    "GALAXY",
    "QSO"
]

vi_std_comments = [
    # Standardized VI comments
    "Broad absorption line quasar (BAL)",
    "Damped Lyman-alpha system (DLA)",
    "Two objects in spectrum",
    "Blazar"
]

_resource_cache = {'templates': None, 'js': None}


def get_resources(filetype):
    """Find all HTML template or JavaScript files in the package.

    Caches the results for quick access.

    Parameters
    ----------
    filetype : {'templates', 'js'}
        The type of file resource needed.

    Returns
    -------
    :class:`dict`
        A dictionary mapping filename to the contents of the file.

    Raises
    ------
    ValueError
        If `filetype` is unknown.
    """
    global _resource_cache
    if filetype not in _resource_cache:
        raise ValueError("Unknown filetype '{0}' for get_resources()!".format(filetype))
    if _resource_cache[filetype] is None:
        _resource_cache[filetype] = dict()
        for f in resource_listdir('prospect', filetype):
            if not f.startswith("."):
                _resource_cache[filetype][f] = resource_string('prospect', filetype + '/' + f).decode('utf-8')
    return _resource_cache[filetype]


def match_catalog_to_spectra(zcat_in, spectra, return_index=False):
    """ Creates a subcatalog, matching a set of DESI spectra

    Parameters
    ----------
    zcat_in : :class:`~astropy.table.Table`, with TARGETID keys
    spectra : :class:`~desispec.spectra.Spectra`
    return_index : :class:`bool`, optional
        If ``True``, returns the list of indices in zcat_in which match spectra

    Returns
    -------
    :class:`~astropy.table.Table`
        A subtable of zcat_in, with rows matching input spectra's TARGETIDs
        If return_index is ``True``, returns (subtable, list of indices)

    Raises
    ------
    RuntimeError
        If a unique row in zcat_in is not found matching each of spectra's TARGETIDs
    """

    if zcat_in is None : return None

    zcat_out = Table(dtype=zcat_in.dtype)
    index_list = list()
    for i_spec in range(spectra.num_spectra()) :
        ww, = np.where((zcat_in['TARGETID'] == spectra.fibermap['TARGETID'][i_spec]))
        if len(ww)<1 :
            raise RuntimeError("No entry in zcat_in for TARGETID "+str(spectra.fibermap['TARGETID'][i_spec]))
        elif len(ww)>1 :
            raise RuntimeError("Several entries in zcat_in for TARGETID "+str(spectra.fibermap['TARGETID'][i_spec]))
        zcat_out.add_row(zcat_in[ww[0]])
        index_list.append(ww[0])
    if return_index:
        return (zcat_out, index_list)
    else:
        return zcat_out


def match_redrockfile_to_spectra(redrockfile, spectra, Nfit=None):
    """ Creates a Table from a Redrock output fit, matching a list of DESI spectra.

    Parameters
    ----------
    redrockfile : :class:`str`, filename for the Redrock output file (.h5 file)
    spectra : :class:`~desispec.spectra.Spectra`
    Nfit : :class:`int`, optional
        Number of best-fits to store in output Table. By default, store all fits available in Redrock file

    Returns
    -------
    :class:`~astropy.table.Table`
        Table with the following columns: TARGETID, CHI2, DELTACHI2, COEFF, Z, ZERR, ZWARN, SPECTYPE, SUBTYPE.
        The rows are matched to spectra's TARGETIDs

    Raises
    ------
    RuntimeError
        If a set of Nfit rows in redrockfile is not found matching each of spectra's TARGETIDs
    """

    dummy, rr_table = redrock.results.read_zscan(redrockfile)
    rr_targets = rr_table['targetid']
    if Nfit is None:
        ww, = np.where( (rr_targets == rr_targets[0]) )
        Nfit = len(ww)
    matched_redrock_cat = Table(
        dtype=[('TARGETID', '<i8'), ('CHI2', '<f8', (Nfit,)),
               ('DELTACHI2', '<f8', (Nfit,)), ('COEFF', '<f8', (Nfit,10,)),
               ('Z', '<f8', (Nfit,)), ('ZERR', '<f8', (Nfit,)),
               ('ZWARN', '<i8', (Nfit,)), ('SPECTYPE', '<U6', (Nfit,)), ('SUBTYPE', '<U2', (Nfit,))])

    for i_spec in range(spectra.num_spectra()):
        ww, = np.where((rr_targets == spectra.fibermap['TARGETID'][i_spec]))
        if len(ww)<Nfit :
            raise RuntimeError("Redrock table cannot match spectra with "+str(Nfit)+" best fits")
        ind = np.argsort(rr_table[ww]['chi2'])[0:Nfit] # Sort fit results by chi2 (independently of spectype)
        sub_table = rr_table[ww][ind]
        the_entry = [ spectra.fibermap['TARGETID'][i_spec] ]
        for redrock_key in ['chi2', 'deltachi2', 'coeff', 'z', 'zerr', 'zwarn', 'spectype', 'subtype']:
            the_entry.append(sub_table[redrock_key])
        matched_redrock_cat.add_row(the_entry)

    return matched_redrock_cat


def create_zcat_from_redrock_cat(redrock_cat, fit_num=0):
    """ Extract a single fit catalog from a redrock catalog containing several fit results per TARGETID

    Parameters
    ----------
    redrock_cat : :class:`~astropy.table.Table`
        Catalog with rows as defined in `match_redrockfile_to_spectra()`
    fit_num : :class:`int`, optional
        The (fit_num)th fit in redrock_cat is extracted (default: 0 ie. redrock's best fit)

    Returns
    -------
    :class:`~astropy.table.Table`
        Table with the following columns: TARGETID, CHI2, COEFF, Z, ZERR, ZWARN, SPECTYPE, SUBTYPE, DELTACHI2.
        ie. identical to the 'ZBEST' table in Redrock's `zbest` files, a single fit per row.
    """

    rr_cat_num_best_fits = redrock_cat['Z'].shape[1]
    if (fit_num >= rr_cat_num_best_fits):
        raise ValueError("fit_num too large wrt redrock_cat")
    zcat_dtype=[('TARGETID', '<i8'), ('CHI2', '<f8'), ('COEFF', '<f8', (10,)),
                ('Z', '<f8'), ('ZERR', '<f8'), ('ZWARN', '<i8'),
                ('SPECTYPE', '<U6'), ('SUBTYPE', '<U2'), ('DELTACHI2', '<f8')]
    zcat_out = Table( data=np.zeros(len(redrock_cat), dtype=zcat_dtype) )
    zcat_out['TARGETID'] = redrock_cat['TARGETID']
    for key in ['CHI2', 'DELTACHI2', 'COEFF', 'SPECTYPE', 'SUBTYPE', 'Z', 'ZERR', 'ZWARN']:
        zcat_out[key] = redrock_cat[key][:,fit_num]

    return zcat_out


def get_subset_label(subset, dirtree_type):
    if dirtree_type=='cumulative':
        label = 'thru'+subset
    elif dirtree_type=='perexp':
        label = 'exp'+subset
    elif dirtree_type=='pernight':
        label = subset
    elif dirtree_type=='exposures':
        label = subset
    else:
        raise ValueError("Unrecognized value for dirtree_type.")
    return label

def create_subsetdb(datadir, dirtree_type=None, spectra_type='coadd', tiles=None, nights=None, expids=None, petals=None, with_zcat=True):
    """Create a 'mini-db' of DESI spectra files, in a given directory tree.

    Supports tile-based and exposure-based directory trees for daily, andes, ... to denali.
    This routine does not open any file, it just checks they exist.

    Parameters
    ----------
    datadir : :class:`string`
        No description provided.
    dirtree_type : :class:`string`
        The directory tree and file names must match the types listed in the notes below.
    spectra_type : :class:`string`, optional
        [c/s]frames are only supported when dirtree_type='exposures'
    petals : :class:`list`, optional
        Filter a set of petal numbers.
    tiles : :class:`list`, optional
        Filter a list of tiles.
    nights : :class:`list`, optional
        Filter a list of nights (only if dirtree_type='pernight' or 'exposures').
    expids : :class:`list`, optional
        Filter a list of exposures (only if dirtree_type='perexp' or 'exposures').
    with_zcat : :class:`bool`, optional
        If True, filter spectra for which a 'zbest' file exists at the same location.

    Returns
    -------
    :class:`dict`
        Content of the 'mini-db':

        - if dirtree_type='exposures': [ {'dataset':night, 'subset':expid, 'petals':[list of petals]}]
        - if dirtree_type='perexp':    [ {'dataset':tile, 'subset':expid, 'petals':[list of petals]}]
        - else:                        [ {'dataset':tile, 'subset':night, 'petals':[list of petals]}]

    Notes
    -----
    * `dirtree_type` must be one of the following:

       - ``dirtree_type='pernight'``: ``{datadir}/{tileid}/{night}/{spectra_type}-{petal}-{tile}-{night}.fits``
       - ``dirtree_type='perexp'``: ``{datadir}/{tileid}/{expid}/{spectra_type}-{petal}-{tile}-exp{expid}.fits``
       - ``dirtree_type='cumulative'``: ``{datadir}/{tileid}/{night}/{spectra_type}-{petal}-{tile}-thru{night}.fits``
       - ``dirtree_type='exposures'``: ``{datadir}/{night}/{expid}/{spectra_type}-{band}{petal}-{expid}.fits``
       - Note that 'perexp' and 'exposures' are different.
       - To use blanc/cascades 'all' (resp 'deep') coadds, use dirtree_type='pernight' and nights=['all'] (resp ['deep']).

    """
    # TODO support (everest) healpix-based directory trees

    if ( (nights is not None and dirtree_type!='pernight' and dirtree_type!='exposures')
        or (expids is not None and dirtree_type!='perexp' and dirtree_type!='exposures') ):
        raise ValueError('Nights/expids option is incompatible with dirtree_type.')
    if dirtree_type == 'exposures':
        if spectra_type not in ['frame', 'cframe', 'sframe']:
            raise ValueError('Unsupported spectra_type: '+spectra_type)
        if with_zcat:
            raise ValueError('Cannot filter zbest files when dirtree_type=exposures')
    else:
        if spectra_type not in ['coadd', 'spectra']:
            raise ValueError('Unsupported spectra_type: '+spectra_type)
    if petals is None:
        petals = [str(i) for i in range(10)]

    #- 'datasets': first level in the explored directory tree
    if dirtree_type == 'exposures':
        datasets = nights
    else:
        datasets = tiles
    if datasets is None:
        datasets = os.listdir(datadir)
    else :
        if not all(x in os.listdir(datadir) for x in datasets):
            raise RuntimeError('Some tile[s]/nights[s] were not found in directory tree.')

    subsetdb = list()
    for dataset in datasets:
        #- 'subsets': second level in the explored directory tree
        all_subsets = os.listdir(os.path.join(datadir,dataset))
        if (nights is not None) and (dirtree_type!='exposures'):
            all_subsets = [ x for x in all_subsets if x in nights ]
        elif expids is not None:
            all_subsets = [ x for x in all_subsets if x in expids ]
        else:
            #- No subset selection, but we discard subdirectories with non-decimal names
            all_subsets = [ x for x in all_subsets if x.isdecimal() ]
        for subset in all_subsets:
            existing_petals = []
            for petal in petals:
                subset_label = get_subset_label(subset, dirtree_type)
                if dirtree_type == 'exposures':
                    spectra_fnames = [ spectra_type+'-'+band+petal+'-'+subset_label+'.fits' for band in ['b', 'r', 'z'] ]
                    if all([os.path.isfile(os.path.join(datadir, dataset, subset, x)) for x in spectra_fnames]):
                        existing_petals.append(petal)
                else:
                    file_label = '-'.join([petal, dataset, subset_label])
                    spectra_fname = os.path.join(datadir, dataset, subset, spectra_type+'-'+file_label+'.fits')
                    zcat_fname = os.path.join(datadir, dataset, subset, 'zbest-'+file_label+'.fits')
                    if os.path.isfile(spectra_fname) and ( (not with_zcat) or os.path.isfile(zcat_fname) ):
                        existing_petals.append(petal)
            if len(existing_petals)>0:
                subsetdb.append( {'dataset':dataset, 'subset':subset, 'petals':existing_petals} )

    return subsetdb

def create_targetdb(datadir, subsetdb, dirtree_type=None):
    """Create a "mini-db" of DESI targetids.

        To do so, `zbest` files are read (faster than reading spectra).

        Parameters
        ----------
        datadir : :class:`string`
            No description provided.
        subsetdb: :class:`list`
            List of spectra subsets, as produced by `create_subsetdb`.
            Format: [ {'tile':tile, 'subset':subset, 'petal':petal} ]
        dirtree_type : :class:`string`
            See documentation in `create_subsetdb`.
            dirtree_type='exposures' is not supported here (no zbest file available in that case).
            Tile-based directory trees for daily, andes, ... to denali are supported.

        Returns
        -------
        :class:`dict`
            Content of the "mini-db": { (tile, subset, petal): [list of TARGETIDs] } where subset is a night or expid.
    """
    if dirtree_type=='exposures':
        raise ValueError("dirtree_type='exposures' is not supported in `create_targetdb`")
    targetdb = dict()
    for the_entry in subsetdb:
        subset_label = get_subset_label(the_entry['subset'], dirtree_type)
        for petal in the_entry['petals']:
            fname = os.path.join(datadir, the_entry['dataset'], the_entry['subset'],
                                 'zbest-'+petal+'-'+the_entry['dataset']+'-'+subset_label+'.fits')
            targetids = np.unique(Table.read(fname,'ZBEST')['TARGETID'])
            targetdb[ (the_entry['dataset'], the_entry['subset'], petal) ] = np.array(targetids, dtype='int64')

    return targetdb


def load_spectra_zcat_from_targets(targetids, datadir, targetdb, dirtree_type='pernight', with_redrock=False, with_redrock_version=True):
    """Get spectra, 'ZBEST' catalog and optional full Redrock catalog matched to a set of DESI TARGETIDs.

    This works using a "mini-db" of targetids, as returned by `create_targetdb()`.
    The outputs of this utility can be used directly by `viewer.plotspectra()`, to inspect a given list of targetids.
    Output spectra/catalog(s) are sorted according to the input target list.
    When several spectra are available for a given TARGETID, they are all included in the output, in random order.

    Parameters
    ----------
    targetids : :class:`list` or :class:`numpy.ndarray`
        List of TARGETIDs, must be int64.
    datadir : :class:`string`
        No description provided.
    dirtree_type : :class:`string`
        The directory tree and file names must match the types listed in the notes below.
    targetdb : :class:`dict`
        Content of the "mini-db": { (tile, subset, petal): [list of TARGETIDs] } where subset is a night or expid.
    with_redrock : :class:`bool`, optional
        If `True`, Redrock output files (.h5 files) are also read
    with_redrock_version : :class:`bool`, optional
        If `True`, a column 'RRVER' is appended to the output redshift catalog, as given by HDU0 in `ZBEST` files.
        This is used by `viewer.plotspectra()` to track Redrock version in visual inspection files.

    Returns
    -------
    :func:`tuple`
        If with_redrock is `False` (default), returns (spectra, zcat), where spectra is `~desispec.spectra.Spectra`
        and zcat is `~astropy.table.Table`.
        If with_redrock is `True`, returns (spectra, zcat, redrockcat) where redrockcat is `~astropy.table.Table`.

    Notes
    -----
    * `dirtree_type` must be one of the following, for "coadd", "zbest" and "redrock" files:

      - ``dirtree_type='pernight'``: ``{datadir}/{tileid}/{night}/zbest-{petal}-{tile}-{night}.fits``
      - ``dirtree_type='perexp'``: ``{datadir}/{tileid}/{expid}/zbest-{petal}-{tile}-exp{expid}.fits``
      - ``dirtree_type='cumulative'``: ``{datadir}/{tileid}/{night}/zbest-{petal}-{tile}-thru{night}.fits``
      - To use blanc/cascades 'all' (resp 'deep') coadds, use dirtree_type='pernight' and nights=['all'] (resp 'deep')

    """

    targetids = np.asarray(targetids)
    if targetids.dtype not in ['int64', 'i8', '>i8']:
        raise TypeError('TARGETIDs should be int64')
    spectra = None
    ztables, rrtables = [], []

    for tile, subset, petal in targetdb.keys():
        targets_subset = set(targetdb[tile, subset, petal])
        targets_subset = targets_subset.intersection(set(targetids))

        # Load spectra for that tile-subset-petal only if one or more target(s) are in the list
        if len(targets_subset)>0 :
            subset_label = get_subset_label(subset, dirtree_type)
            file_label = '-'.join([petal, tile, subset_label])
            the_path = os.path.join(datadir, tile, subset)
            the_spec = desispec.io.read_spectra(os.path.join(the_path, "coadd-"+file_label+".fits"))
            the_spec = the_spec.select(targets=sorted(targets_subset))
            the_zcat = Table.read(os.path.join(the_path, "zbest-"+file_label+".fits"), 'ZBEST')
            if with_redrock_version:
                hdulist = astropy.io.fits.open(os.path.join(the_path, "zbest-"+file_label+".fits"))
                the_zcat['RRVER'] = hdulist[hdulist.index_of('PRIMARY')].header['RRVER']
            the_zcat = match_catalog_to_spectra(the_zcat, the_spec)
            ztables.append(the_zcat)
            if with_redrock:
                rrfile = os.path.join(the_path, "redrock-"+file_label+".h5")
                the_rrcat = match_redrockfile_to_spectra(rrfile, the_spec, Nfit=None)
                rrtables.append(the_rrcat)
            if spectra is None:
                spectra = the_spec
            else:
                #- Still use update() instead of stack(), to handle case when fibermaps differ in different files.
                spectra.update(the_spec)

    #- Sort according to input target list. Check if all targets were found in spectra
    tids_spectra = spectra.fibermap['TARGETID']
    sorted_indices = []
    for target in targetids:
        w, = np.where(tids_spectra == target)
        sorted_indices.extend(w)
        if len(w)==0:
            print("Warning! TARGETID not found:", target)
    assert(len(tids_spectra)==len(sorted_indices)) # check, should always be true
    spectra = spectra[ sorted_indices ]

    zcat = vstack(ztables)
    zcat = zcat[ sorted_indices ]
    if with_redrock:
        rrcat = vstack(rrtables)
        rrcat = rrcat[ sorted_indices ]
        return (spectra, zcat, rrcat)
    else:
        return (spectra, zcat)


def frames2spectra(frames, nspec=None, startspec=None, with_scores=False, with_resolution_data=False):
    """Convert list of frames into DESI Spectra object

    Parameters
    ----------
    frames : :class:`list`
        A list of :class:`~desispec.frame.Frame`.
    nspec : :class:`int`, optional
        No description provided.
    startspec : :class:`int`, optional
        If nspec is set, only spectra in range [startspec:nspec+startspec] are kept
    with_scores : :class:`bool`, optional
        If `True`, include merged scores from input frames
    with_resolution_data : :class:`bool`, optional
        If `True`, include frames.resolution_data

    Returns
    -------
    :class:`~desispec.spectra.Spectra`
        No description provided.
    """

    bands = list()
    wave = dict()
    flux = dict()
    ivar = dict()
    mask = dict()
    res = dict()

    for fr in frames:
        fibermap = fr.fibermap
        band = fr.meta['CAMERA'][0]
        bands.append(band)
        wave[band] = fr.wave
        flux[band] = fr.flux
        ivar[band] = fr.ivar
        mask[band] = fr.mask
        res[band] = fr.resolution_data
        if nspec is not None :
            if startspec is None : startspec = 0
            flux[band] = flux[band][startspec:nspec+startspec]
            ivar[band] = ivar[band][startspec:nspec+startspec]
            mask[band] = mask[band][startspec:nspec+startspec]
            res[band] = res[band][startspec:nspec+startspec,:,:]
            fibermap = fr.fibermap[startspec:nspec+startspec]

    merged_scores = None
    if with_scores :
        scores_columns = frames[0].scores.columns
        for i in range(1,len(frames)) :
            scores_columns += frames[i].scores.columns
        merged_scores = astropy.io.fits.FITS_rec.from_columns(scores_columns)

    if not with_resolution_data : res = None

    spectra = desispec.spectra.Spectra(
        bands, wave, flux, ivar, mask, fibermap=fibermap, meta=fr.meta, scores=merged_scores, resolution_data=res
    )
    return spectra


def metadata_selection(spectra, mask=None, mask_type=None, gmag_range=None, rmag_range=None, chi2_range=None, snr_range=None, clean_fiberstatus=False, with_dirty_mask_merge=False, zcat=None, log=None):
    """Simple selection of DESI spectra based on various metadata.

    Filtering based on the logical AND of requested selection criteria.
    Note: use X_range=[min, None] to filter X > min, X_range=[None, max] to filter X < max

    Parameters
    ----------
    spectra : :class:`~desispec.spectra.Spectra`
        No description provided.
    mask : :class:`string`, optional
        DESI targeting mask to select, eg 'ELG'. Requires to set mask_type.
    mask_type : :class:`string`, optional
        DESI targeting mask category, currently supported: 'DESI_TARGET', 'BGS_TARGET',
        'MWS_TARGET', 'SECONDARY_TARGET', 'CMX_TARGET', 'SV[1/2/3]_DESI_TARGET', 'SV[1/2/3]_BGS_TARGET',
        'SV[1/2/3]_MWS_TARGET', 'SV[1/2/3]_SCND_TARGET'.
    with_dirty_mask_merge : :class:`bool`, optional
        Option for specific targeting mask selection in early CMX data, see code...
    gmag_range : :class:`list`
        g magnitude range to select, gmag_range = [gmag_min, gmag_max]
    rmag_range : :class:`list`
        r magnitude range to select, rmag_range = [rmag_min, rmag_max]
    snr_range : :class:`list`
        SNR range to select, snr_range = [snr_min, snr_max].
        This filter applies on all B, R and Z bands, from scores.MEDIAN_COADD_SNR_band, or
        scores.MEDIAN_CALIB_SNR_band if the former is not found.
    chi2_range : :class:`list`
        chi2 range to select, chi2_range = [chi2_min, chi2_max]. Requires to set zcat.
    clean_fiberstatus : :class:`bool`
        if True, remove spectra with FIBERSTATUS!=0
    zcat : :class:`~astropy.table.Table`
        catalog with chi2 information, must be matched to spectra (needed for chi2_range filter).
    log : optional log.


    Returns
    -------
    :class:`~desispec.spectra.Spectra`
        No description provided.
    """
    keep = np.ones(len(spectra.fibermap), bool)

    #- SNR selection
    if (snr_range is not None) and (snr_range!=[None, None]):
        #- If a bound is set to None, replace by +-np.inf
        if snr_range[0]==None:
            snr_range[0] = -np.inf
        if snr_range[1]==None:
            snr_range[1] = np.inf
        if len(snr_range)!=2 or snr_range[1]<snr_range[0]:
            raise ValueError("Wrong input snr_range")
        if spectra.scores is None:
            raise RuntimeError('No scores in spectra: cannot select on SNR')
        snr_var = 'MEDIAN_COADD_SNR'
        if snr_var+'_B' not in spectra.scores.keys():
            snr_var = 'MEDIAN_CALIB_SNR'
        for band in ['B','R','Z'] :
            keep_snr = ( (spectra.scores[snr_var+'_'+band]>snr_range[0]) &
                       (spectra.scores[snr_var+'_'+band]<snr_range[1]) )
            if np.all(~keep_snr):
                if log is not None :
                    log.info(" * No spectra with MEDIAN_CALIB_SNR_"+band+" in requested range")
                return None
            else :
                keep = ( keep & keep_snr )

    #- Target mask selection
    if mask is not None :
        if not _desitarget_imported:
            raise RuntimeError('desitarget not imported: cannot select on targeting mask')
        if mask_type not in spectra.fibermap.keys():
            mask_candidates = [x for x in spectra.fibermap.keys() if '_TARGET' in x]
            raise ValueError(mask_type+" is not in spectra.fibermap.\n Hints of available masks: "+(' '.join(mask_candidates)))
        mask_used = supported_desitarget_masks[mask_type]
        if mask not in mask_used.names():
            raise ValueError("requested mask "+mask+" does not match mask_type "+mask_type)
        keep_mask = (spectra.fibermap[mask_type] & mask_used[mask]) > 0   # boolean array
        if mask_type == 'CMX_TARGET' and with_dirty_mask_merge:
            #- Self-explanatory... only for fast VI of minisv
            mask2 = None
            if mask in ['SV0_QSO', 'SV0_ELG', 'SV0_LRG']: mask2 = mask.replace('SV0','MINI_SV')
            if mask == 'SV0_BGS': mask2 = 'MINI_SV_BGS_BRIGHT'
            if mask in ['SV0_STD_FAINT', 'SV0_STD_BRIGHT']: mask2 = mask.replace('SV0_','')
            if mask2 is not None:
                keep_mask = ( (spectra.fibermap[mask_type] & mask_used[mask]) |
                            (spectra.fibermap[mask_type] & mask_used[mask2]) ) > 0
        if np.all(~keep_mask):
            if log is not None : log.info(" * No spectra with mask "+mask)
            return None
        else :
            keep = ( keep & keep_mask )

    #- Photometry selection
    if (gmag_range is not None) and (gmag_range!=[None, None]):
        if gmag_range[0]==None:
            gmag_range[0] = -np.inf
        if gmag_range[1]==None:
            gmag_range[1] = np.inf
        if len(gmag_range)!=2 or gmag_range[1]<gmag_range[0]:
            raise ValueError("Wrong input gmag_range")
        gmag = np.zeros(spectra.num_spectra())
        w, = np.where( (spectra.fibermap['FLUX_G']>0) )
        gmag[w] = -2.5*np.log10(spectra.fibermap['FLUX_G'][w])+22.5
        if 'MW_TRANSMISSION_G' in spectra.fibermap.keys():
            w, = np.where( (spectra.fibermap['FLUX_G']>0) & (spectra.fibermap['MW_TRANSMISSION_G']>0) )
            gmag[w] = -2.5*np.log10(spectra.fibermap['FLUX_G'][w]/spectra.fibermap['MW_TRANSMISSION_G'][w])+22.5
        keep_gmag = ( (gmag>gmag_range[0]) & (gmag<gmag_range[1]) )
        if np.all(~keep_gmag):
            if log is not None : log.info(" * No spectra with g_mag in requested range")
            return None
        else :
            keep = ( keep & keep_gmag )

    if (rmag_range is not None) and (rmag_range!=[None, None]):
        if rmag_range[0]==None:
            rmag_range[0] = -np.inf
        if rmag_range[1]==None:
            rmag_range[1] = np.inf
        if len(rmag_range)!=2 or rmag_range[1]<rmag_range[0]:
            raise ValueError("Wrong input rmag_range")
        rmag = np.zeros(spectra.num_spectra())
        w, = np.where( (spectra.fibermap['FLUX_R']>0) )
        rmag[w] = -2.5*np.log10(spectra.fibermap['FLUX_R'][w])+22.5
        if 'MW_TRANSMISSION_R' in spectra.fibermap.keys():
            w, = np.where( (spectra.fibermap['FLUX_R']>0) & (spectra.fibermap['MW_TRANSMISSION_R']>0) )
            rmag[w] = -2.5*np.log10(spectra.fibermap['FLUX_R'][w]/spectra.fibermap['MW_TRANSMISSION_R'][w])+22.5
        keep_rmag = ( (rmag>rmag_range[0]) & (rmag<rmag_range[1]) )
        if np.all(~keep_rmag):
            if log is not None : log.info(" * No spectra with r_mag in requested range")
            return None
        else :
            keep = ( keep & keep_rmag )

    #- Chi2 selection
    if (chi2_range is not None) and (chi2_range!=[None, None]):
        if chi2_range[0]==None:
            chi2_range[0] = -np.inf
        if chi2_range[1]==None:
            chi2_range[1] = np.inf
        if len(chi2_range)!=2 or chi2_range[1]<chi2_range[0]:
            raise ValueError("Wrong input chi2_range")
        if np.any(zcat['TARGETID'] != spectra.fibermap['TARGETID']) :
            raise RuntimeError('zcat and spectra do not match (different targetids)')
        keep_chi2 = ( (zcat['DELTACHI2']>chi2_range[0]) & (zcat['DELTACHI2']<chi2_range[1]) )
        if np.all(~keep_chi2):
            if log is not None : log.info(" * No target in this pixel with DeltaChi2 in requested range")
            return None
        else :
            keep = ( keep & keep_chi2 )

    #- Fiberstatus selection
    if clean_fiberstatus and ('FIBERSTATUS' in spectra.fibermap.keys()):
        keep = ( keep & (spectra.fibermap['FIBERSTATUS']==0) )

    return spectra[keep]


def _coadd(wave, flux, ivar, rdat):
    '''Return weighted coadd of spectra

    Parameters
    ----------
    wave : array-like
        1D[nwave] array of wavelengths.
    flux : array-like
        2D[nspec, nwave] array of flux densities.
    ivar : array-like
        2D[nspec, nwave] array of inverse variances of `flux`.
    rdat : array-like
        3D[nspec, ndiag, nwave] sparse diagonals of resolution matrix.

    Returns
    -------
    :class:`tuple`
        The coadded spectrum (wave, outflux, outivar, outrdat).
    '''
    nspec, nwave = flux.shape
    unweightedflux = np.zeros(nwave, dtype=flux.dtype)
    weightedflux = np.zeros(nwave, dtype=flux.dtype)
    weights = np.zeros(nwave, dtype=flux.dtype)
    outrdat = np.zeros(rdat[0].shape, dtype=rdat.dtype)
    for i in range(nspec):
        unweightedflux += flux[i]
        weightedflux += flux[i] * ivar[i]
        weights += ivar[i]
        outrdat += rdat[i] * ivar[i]

    isbad = (weights == 0)
    outflux = weightedflux / (weights + isbad)
    outflux[isbad] = unweightedflux[isbad] / nspec

    outrdat /= (weights + isbad)
    outivar = weights

    return wave, outflux, outivar, outrdat

def coadd_targets(spectra, targetids=None):
    '''
    Coadds individual exposures of the same targets; returns new Spectra object

    Parameters
    ----------
    spectra : :class:`desispec.spectra.Spectra`
    targetids : array-like, optional
        Subset of target IDs to keep.

    Returns
    -------
    :class:`desispec.spectra.Spectra`
        Where individual spectra of each target have been combined into a
        single spectrumper camera.

    Notes
    -----
    Coadds per camera but not across cameras.
    '''
    if targetids is None:
        targetids = spectra.target_ids()

    #- Create output arrays to fill
    ntargets = spectra.num_targets()
    wave = dict()
    flux = dict()
    ivar = dict()
    rdat = dict()
    if spectra.mask is None:
        mask = None
    else:
        mask = dict()
    for channel in spectra.bands:
        wave[channel] = spectra.wave[channel].copy()
        nwave = len(wave[channel])
        flux[channel] = np.zeros((ntargets, nwave))
        ivar[channel] = np.zeros((ntargets, nwave))
        ndiag = spectra.resolution_data[channel].shape[1]
        rdat[channel] = np.zeros((ntargets, ndiag, nwave))
        if mask is not None:
            mask[channel] = np.zeros((ntargets, nwave), dtype=spectra.mask[channel].dtype)

    #- Loop over targets, coadding all spectra for each target
    fibermap = Table(dtype=spectra.fibermap.dtype)
    for i, targetid in enumerate(targetids):
        ii = np.where(spectra.fibermap['TARGETID'] == targetid)[0]
        fibermap.add_row(spectra.fibermap[ii[0]])
        for channel in spectra.bands:
            if len(ii) > 1:
                outwave, outflux, outivar, outrdat = _coadd(
                    spectra.wave[channel],
                    spectra.flux[channel][ii],
                    spectra.ivar[channel][ii],
                    spectra.resolution_data[channel][ii]
                    )
                if mask is not None:
                    outmask = spectra.mask[channel][ii[0]]
                    for j in range(1, len(ii)):
                        outmask |= spectra.mask[channel][ii[j]]
            else:
                outwave, outflux, outivar, outrdat = (
                    spectra.wave[channel],
                    spectra.flux[channel][ii[0]],
                    spectra.ivar[channel][ii[0]],
                    spectra.resolution_data[channel][ii[0]]
                    )
                if mask is not None:
                    outmask = spectra.mask[channel][ii[0]]

            flux[channel][i] = outflux
            ivar[channel][i] = outivar
            rdat[channel][i] = outrdat
            if mask is not None:
                mask[channel][i] = outmask

    return desispec.spectra.Spectra(spectra.bands, wave, flux, ivar,
            mask=mask, resolution_data=rdat, fibermap=fibermap,
            meta=spectra.meta)

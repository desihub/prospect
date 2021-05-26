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
    supported_desitarget_masks = {
        'DESI_TARGET': desi_mask,
        'BGS_TARGET': bgs_mask,
        'MWS_TARGET': mws_mask,
        'SECONDARY_TARGET': scnd_mask, # To confirm !
        'CMX_TARGET': cmx_mask,
        'SV1_DESI_TARGET': sv1_desi_mask,
        'SV1_BGS_TARGET': sv1_bgs_mask,
        'SV1_MWS_TARGET': sv1_mws_mask,
        'SV1_SCND_TARGET': sv1_scnd_mask,
        }
except ImportError:
    _desitarget_imported = False
    supported_desitarget_masks = dict()

_redrock_imported = True
try:
    import redrock.results
except ImportError:
    _redrock_imported = False


# from prospect import mycoaddcam  # Does not appear to be used in this module.
from prospect import myspecselect, myspecupdate

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


def match_zcat_to_spectra(zcat_in, spectra) :
    '''
    zcat_in : astropy Table from redshift fitter
    - creates a new astropy Table whose rows match the targetids of input spectra
    - also returns the corresponding list of indices
    - for each targetid, a unique row in zcat_in must exist.
    TODO : maybe rename this fct ? match_table_to_spectra ?
    => it also works whatever kind of input zcat : just has to be a table with 'TARGETID' key
    => in particular it's useful for "redrock_cat" tables
    '''

    if zcat_in is None : return None

    zcat_out = Table(dtype=zcat_in.dtype)
    index_list = list()
    for i_spec in range(spectra.num_spectra()) :
        ww, = np.where((zcat_in['TARGETID'] == spectra.fibermap['TARGETID'][i_spec]))
        if len(ww)<1 :
                raise RuntimeError("No zcat entry for target "+str(spectra.fibermap['TARGETID'][i_spec]))
        if len(ww)>1 :
            raise RuntimeError("Several zcat entries for target "+str(spectra.fibermap['TARGETID'][i_spec]))
        zcat_out.add_row(zcat_in[ww[0]])
        index_list.append(ww[0])
    return (zcat_out, index_list)


def match_redrock_zfit_to_spectra(redrockfile, spectra, Nfit=None) :
    '''
    Read Redrock file, and return astropy Table of best fits matched to the targetids of input spectra
    - for each target, store arrays chi2[Nfit], coeff[Nfit], z[Nfit], spectype[Nfit], subtype[Nfit]
    - if Nfit is None: take all available fits
    '''

    dummy, rr_table = redrock.results.read_zscan(redrockfile)
    rr_targets = rr_table['targetid']
    if Nfit is None :
        ww, = np.where( (rr_targets == rr_targets[0]) )
        Nfit = len(ww)
    matched_redrock_cat = Table(dtype=[('TARGETID', '<i8'), ('CHI2', '<f8', (Nfit,)), ('DELTACHI2', '<f8', (Nfit,)), ('COEFF', '<f8', (Nfit,10,)), ('Z', '<f8', (Nfit,)), ('ZERR', '<f8', (Nfit,)), ('ZWARN', '<i8', (Nfit,)), ('SPECTYPE', '<U6', (Nfit,)), ('SUBTYPE', '<U2', (Nfit,))])

    for i_spec in range(spectra.num_spectra()) :
        ww, = np.where((rr_targets == spectra.fibermap['TARGETID'][i_spec]))
        if len(ww)<Nfit :
            raise RuntimeError("redrock table cannot match spectra with "+str(Nfit)+" best fits")
        ind = np.argsort(rr_table[ww]['chi2'])[0:Nfit]
        sub_table = rr_table[ww][ind]
        the_entry = [ spectra.fibermap['TARGETID'][i_spec] ]
        the_entry.append(sub_table['chi2'])
        the_entry.append(sub_table['deltachi2'])
        the_entry.append(sub_table['coeff'])
        the_entry.append(sub_table['z'])
        the_entry.append(sub_table['zerr'])
        the_entry.append(sub_table['zwarn'])
        the_entry.append(sub_table['spectype'])
        the_entry.append(sub_table['subtype'])
        matched_redrock_cat.add_row(the_entry)

    return matched_redrock_cat


def create_zcat_from_redrock_cat(redrock_cat, fit_num=0) :
    '''
    TODO change name zcat -> zbest_cat ?
    Extract a z catalog from redrock catalog produced in match_redrock_zfit_to_spectra()
    The z catalog has one fit per targetid, corresponding to the (fit_num)th best fit
    '''

    rr_cat_num_best_fits = redrock_cat['Z'].shape[1]
    if (fit_num >= rr_cat_num_best_fits) : raise ValueError("fit_num too large wrt redrock_cat")
    zcat_dtype=[('TARGETID', '<i8'), ('CHI2', '<f8'), ('COEFF', '<f8', (10,)), ('Z', '<f8'), ('ZERR', '<f8'), ('ZWARN', '<i8'), ('SPECTYPE', '<U6'), ('SUBTYPE', '<U2'), ('DELTACHI2', '<f8')]
    zcat_out = Table( data=np.zeros(len(redrock_cat), dtype=zcat_dtype) )
    zcat_out['TARGETID'] = redrock_cat['TARGETID']
    for key in [ 'CHI2', 'DELTACHI2', 'COEFF', 'SPECTYPE', 'SUBTYPE', 'Z', 'ZERR', 'ZWARN'] :
        zcat_out[key] = redrock_cat[key][:,fit_num]

    return zcat_out

def make_targetdict(tiledir, petals=[str(i) for i in range(10)], tiles=None, nights=None, cumulative=False) :
    '''
    Small homemade/hack utility (based on Anand's hack)
    Makes "mini-db" of targetids. It basically reads zbest, so it's reasonably fast
    Adapted to the directory structure tiledir/{tile}/{night}/{cframe/zbest/...}
    - tiles: optional list of tiles (all of them must exist in tiledir)
    - nights: optional list of nights
    '''

    if ( ('cumulative' in tiledir) and not cumulative) or ( ('cumulative' not in tiledir) and cumulative):
        print('Warning: cumulative option does not seem to fit tiledir=',tiledir)

    target_dict = {}
    if tiles is None :
        tiles = os.listdir(tiledir)
    else :
        alltiles = os.listdir(tiledir)
        assert all(x in alltiles for x in tiles)

    for tile in tiles :
        the_nights = os.listdir(os.path.join(tiledir,tile))
        if nights is not None :
            the_nights = [ x for x in the_nights if x in nights ]
        for night in the_nights :
            target_dict[tile+"-"+night] = {}
            # exposures (included in dict only if 'cframe-b' files are present - which is not the case for "deep" coadds in blanc)
            fns = np.sort(glob.glob(os.path.join(tiledir,tile,night,'cframe-b'+petals[0]+'-????????.fits')))
            if len(fns)>0 :
                target_dict[tile+"-"+night]['exps'] = np.array([fn.replace('.','-').split('-')[-2] for fn in fns])
            # targetid, fibres
            targetid,fiber,petal_list = [],[],[]
            for petal in petals:
                night_label = 'thru'+night if cumulative else night
                pp = glob.glob(os.path.join(tiledir,tile,night,'zbest-'+petal+'-'+tile+'-'+night_label+'.fits'))
                if len(pp)>0 :
                    fn        = pp[0]
                    fm = Table.read(fn, 'FIBERMAP')
                    tid, keep = np.unique(fm['TARGETID'], return_index=True)
                    data = fm[keep]
                    #data      = astropy.io.fits.open(fn)['fibermap'].data # Not ok with andes
                    targetid += data['TARGETID'].tolist()
                    fiber    += data['FIBER'].tolist()
                    petal_list    += [petal for i in range(len(data['TARGETID']))]
            target_dict[tile+"-"+night]['targetid'] = np.array(targetid, dtype='int64')
            target_dict[tile+"-"+night]['fiber']    = np.array(fiber)
            target_dict[tile+"-"+night]['petal']    = np.array(petal_list)

    return target_dict

def load_spectra_zcat_from_targets(targets, tiledir, obs_db, with_redrock=False, with_redrock_version=True, cumulative=False) :
    '''
    Creates (spectra,zcat,[redrock_cat]) = (Spectra object, zcatalog Table, [redrock Table]) from a list of targetids
    - targets must be a list of int64
    - obs_db: "mini-db" produced by make_targetdict()
    - with_redrock: if True, also get redrock Table
    - with_redrock_version: if True, add a Column to zcat with RRVER from hdu0 in zbest files
    '''

    if ( ('cumulative' in tiledir) and not cumulative) or ( ('cumulative' not in tiledir) and cumulative):
        print('Warning: cumulative option does not seem to fit tiledir=',tiledir)

    targets = np.asarray(targets)
    if targets.dtype not in ['int64', 'i8', '>i8'] :
        raise TypeError('Targetids should be int64.')
    spectra = None
    ztables, rrtables = [], []

    for tile_night in obs_db.keys() :
        petals = np.unique(obs_db[tile_night]['petal'])
        for petal in petals :
            targets_subset = []
            for target in targets :
                w, = np.where( (obs_db[tile_night]['targetid']==target) & (obs_db[tile_night]['petal']==petal) )
                if len(w) == 0 : continue
                if len(w) > 1 :
                    print("Warning ! Several entries in tile/night "+tile_night+" available for target "+str(target))
                targets_subset.append(target)
            # Load spectra for that tile/night/petal only if there's a target from the list
            if len(targets_subset)>0 :
                tile, night = tile_night.split('-')
                the_path = os.path.join(tiledir, tile, night)
                night_label = 'thru'+night if cumulative else night
                file_label = '-'.join([petal, tile, night_label])
                the_spec = desispec.io.read_spectra(os.path.join(the_path,"coadd-"+file_label+".fits"))
                the_spec = myspecselect.myspecselect(the_spec, targets=targets_subset, remove_scores=True)
                the_zcat = Table.read(os.path.join(the_path,"zbest-"+file_label+".fits"),'ZBEST')
                if with_redrock_version:
                    hdulist = astropy.io.fits.open(os.path.join(the_path,"zbest-"+file_label+".fits"))
                    the_zcat['RRVER'] = hdulist[hdulist.index_of('PRIMARY')].header['RRVER']
                the_zcat, dummy = match_zcat_to_spectra(the_zcat, the_spec)
                ztables.append(the_zcat)
                if with_redrock :
                    rrfile = os.path.join(the_path,"redrock-"+file_label+".h5")
                    the_rrcat = match_redrock_zfit_to_spectra(rrfile, the_spec, Nfit=None)
                    rrtables.append(the_rrcat)
                if spectra is None : spectra = the_spec
                else : spectra = myspecupdate.myspecupdate(spectra,the_spec)

    # Check if all targets were found in spectra
    tids_spectra = spectra.fibermap['TARGETID']
    for target in targets :
        if target not in tids_spectra : print("Warning! targetid not found: "+str(target))

    zcat = vstack(ztables)
    if with_redrock :
        rrcat = vstack(rrtables)
        return (spectra, zcat, rrcat)
    else :
        return (spectra,zcat)


def get_y_minmax(pmin, pmax, data, ispec) :
    '''
    Utility, from plotframe
    '''
    dx = np.sort(data[np.isfinite(data)])
    if len(dx)==0 : return (0,0)
    imin = int(np.floor(pmin*len(dx)))
    imax = int(np.floor(pmax*len(dx)))
    if (imax >= len(dx)) : imax = len(dx)-1
    return (dx[imin],dx[imax])



def frames2spectra(frames, nspec=None, startspec=None, with_scores=False, with_resolution_data=False):
    '''Convert input list of Frames into Spectra object
    with_score : if true, propagate scores
    with_resolution_data: if true, propagate resolution
    '''
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


def specviewer_selection(spectra, log=None, mask=None, mask_type=None, gmag_cut=None, rmag_cut=None, chi2cut=None, zbest=None, snr_cut=None, with_dirty_mask_merge=False, remove_scores=False) :
    '''
    Simple sub-selection on spectra based on meta-data.
        Implemented cuts based on : target mask ; photo mag (g, r) ; chi2 from fit ; SNR (in spectra.scores, BRZ)
        - if chi2cut : a catalog zbest must be provided, with entries matching exactly those of spectra
    '''

    # SNR selection
    if snr_cut is not None :
        assert ( (len(snr_cut)==2) and (spectra.scores is not None) )
        for band in ['B','R','Z'] :
            w, = np.where( (spectra.scores['MEDIAN_CALIB_SNR_'+band]>snr_cut[0]) & (spectra.scores['MEDIAN_CALIB_SNR_'+band]<snr_cut[1]) )
            if len(w) == 0 :
                if log is not None : log.info(" * No spectra with MEDIAN_CALIB_SNR_"+band+" in requested range")
                return 0
            else :
                targetids = spectra.fibermap['TARGETID'][w]
                spectra = myspecselect.myspecselect(spectra, targets=targetids, remove_scores=remove_scores)

    # Target mask selection
    if mask is not None :
        assert _desitarget_imported
        if mask_type not in spectra.fibermap.keys():
            raise ValueError("mask_type is not in spectra.fibermap: "+mask_type)
        mask_used = supported_desitarget_masks[mask_type]
        assert ( mask in mask_used.names() )
        w, = np.where( (spectra.fibermap[mask_type] & mask_used[mask]) )
        if mask_type == 'CMX_TARGET' and with_dirty_mask_merge: # Self-explanatory... only for fast VI of minisv
            mask2 = None
            if mask in ['SV0_QSO', 'SV0_ELG', 'SV0_LRG']: mask2 = mask.replace('SV0','MINI_SV')
            if mask == 'SV0_BGS': mask2 = 'MINI_SV_BGS_BRIGHT'
            if mask in ['SV0_STD_FAINT', 'SV0_STD_BRIGHT']: mask2 = mask.replace('SV0_','')
            if mask2 is not None:
                w, = np.where( (spectra.fibermap[mask_type] & mask_used[mask]) |
                             (spectra.fibermap[mask_type] & mask_used[mask2]) )
        if len(w) == 0 :
            if log is not None : log.info(" * No spectra with mask "+mask)
            return 0
        else :
            targetids = spectra.fibermap['TARGETID'][w]
            spectra = myspecselect.myspecselect(spectra, targets=targetids, remove_scores=remove_scores)

    # Photometry selection
    if gmag_cut is not None :
        assert len(gmag_cut)==2 # Require range [gmin, gmax]
        gmag = np.zeros(spectra.num_spectra())
        w, = np.where( (spectra.fibermap['FLUX_G']>0) & (spectra.fibermap['MW_TRANSMISSION_G']>0) )
        gmag[w] = -2.5*np.log10(spectra.fibermap['FLUX_G'][w]/spectra.fibermap['MW_TRANSMISSION_G'][w])+22.5
        w, = np.where( (gmag>gmag_cut[0]) & (gmag<gmag_cut[1]) )
        if len(w) == 0 :
            if log is not None : log.info(" * No spectra with g_mag in requested range")
            return 0
        else :
            targetids = spectra.fibermap['TARGETID'][w]
            spectra = myspecselect.myspecselect(spectra, targets=targetids)
    if rmag_cut is not None :
        assert len(rmag_cut)==2 # Require range [rmin, rmax]
        rmag = np.zeros(spectra.num_spectra())
        w, = np.where( (spectra.fibermap['FLUX_R']>0) & (spectra.fibermap['MW_TRANSMISSION_R']>0) )
        rmag[w] = -2.5*np.log10(spectra.fibermap['FLUX_R'][w]/spectra.fibermap['MW_TRANSMISSION_R'][w])+22.5
        w, = np.where( (rmag>rmag_cut[0]) & (rmag<rmag_cut[1]) )
        if len(w) == 0 :
            if log is not None : log.info(" * No spectra with r_mag in requested range")
            return 0
        else :
            targetids = spectra.fibermap['TARGETID'][w]
            spectra = myspecselect.myspecselect(spectra, targets=targetids, remove_scores=remove_scores)

    # Chi2 selection
    if chi2cut is not None :
        assert len(chi2cut)==2 # Require range [chi2min, chi2max]
        if np.any(zbest['TARGETID'] != spectra.fibermap['TARGETID']) :
            raise RuntimeError('specviewer_selection : zbest and spectra do not match (different targetids)')

        w, = np.where( (zbest['DELTACHI2']>chi2cut[0]) & (zbest['DELTACHI2']<chi2cut[1]) )
        if len(w) == 0 :
            if log is not None : log.info(" * No target in this pixel with DeltaChi2 in requested range")
            return 0
        else :
            targetids = spectra.fibermap['TARGETID'][w]
            spectra = myspecselect.myspecselect(spectra, targets=targetids, remove_scores=remove_scores)

    return spectra


def _coadd(wave, flux, ivar, rdat):
    '''
    Return weighted coadd of spectra

    Parameters
    ----------
    wave : 1D[nwave] array of wavelengths
    flux : 2D[nspec, nwave] array of flux densities
    ivar : 2D[nspec, nwave] array of inverse variances of `flux`
    rdat : 3D[nspec, ndiag, nwave] sparse diagonals of resolution matrix

    Returns
    -------
        coadded spectrum (wave, outflux, outivar, outrdat)
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
    targetids : (optional) array-like subset of target IDs to keep

    Returns
    -------
    coadded_spectra : :class:`desispec.spectra.Spectra` where individual
        spectra of each target have been combined into a single spectrum
        per camera.

    Note: coadds per camera but not across cameras.
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

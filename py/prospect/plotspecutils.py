# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
======================
prospect.plotspecutils
======================

Plot spectra in the form of specutils-compatible objects.
See :mod:`prospect.specutils` for objects and IO routines.
"""
import os, sys
import argparse

import numpy as np
import scipy.ndimage.filters

from astropy.table import Table
import astropy.io.fits
import astropy.convolution
from specutils import Spectrum1D, SpectrumList

import bokeh.plotting as bk
from bokeh.models import ColumnDataSource, CDSView, IndexFilter
from bokeh.models import CustomJS, LabelSet, Label, Span, Legend, Panel, Tabs, BoxAnnotation
from bokeh.models.widgets import (
    Slider, Button, Div, CheckboxGroup, CheckboxButtonGroup, RadioButtonGroup,
    TextInput, Select, DataTable, TableColumn)
import bokeh.layouts as bl
import bokeh.events

_desitarget_imported = True
try:
    from desitarget.targetmask import desi_mask
    from desitarget.cmx.cmx_targetmask import cmx_mask
    from desitarget.sv1.sv1_targetmask import desi_mask as sv1_desi_mask
except ImportError:
    _desitarget_imported = False

from .utilities import get_resources, create_zcat_from_redrock_cat, vi_flags, vi_file_fields, vi_spectypes, vi_std_comments
from .mycoaddcam import coaddcam_prospect
from .plotframes import add_lines, _airtovac

def create_model(spectra, zbest, archetype_fit=False, archetypes_dir=None, template_dir=None):
    '''
    Returns model_wave[nwave], model_flux[nspec, nwave], row matched to zbest,
    which can be in a different order than spectra.
    - zbest must be entry-matched to spectra.
    '''
    import redrock.templates
    from desispec.interpolation import resample_flux

    if archetype_fit:
        from redrock.archetypes import All_archetypes

    if np.any(zbest['TARGETID'] != spectra.fibermap['TARGETID']) :
        raise ValueError('zcatalog and spectra do not match (different targetids)')

    #- Load redrock templates; redirect stdout because redrock is chatty
    saved_stdout = sys.stdout
    sys.stdout = open('/dev/null', 'w')
    try:
        templates = dict()
        for filename in redrock.templates.find_templates(template_dir=template_dir):
            tx = redrock.templates.Template(filename)
            templates[(tx.template_type, tx.sub_type)] = tx
    except Exception as err:
        sys.stdout = saved_stdout
        raise(err)

    sys.stdout = saved_stdout

    #- Empty model flux arrays per band to fill
    model_flux = dict()
    for band in spectra.bands:
        model_flux[band] = np.zeros(spectra.flux[band].shape)

    for i in range(len(zbest)):
        zb = zbest[i]

        if archetype_fit:
            archetypes = All_archetypes(archetypes_dir=archetypes_dir).archetypes
            archetype  = archetypes[zb['SPECTYPE']]
            coeff      = zb['COEFF']

            for band in spectra.bands:
                wave                = spectra.wave[band]
                wavehash            = hash((len(wave), wave[0], wave[1], wave[-2], wave[-1], spectra.R[band].data.shape[0]))
                dwave               = {wavehash: wave}
                mx                  = archetype.eval(zb['SUBTYPE'], dwave, coeff, wave, zb['Z'])
                model_flux[band][i] = spectra.R[band][i].dot(mx)

        else:
            tx    = templates[(zb['SPECTYPE'], zb['SUBTYPE'])]
            coeff = zb['COEFF'][0:tx.nbasis]
            model = tx.flux.T.dot(coeff).T

            for band in spectra.bands:
                mx                  = resample_flux(spectra.wave[band], tx.wave*(1+zb['Z']), model)
                model_flux[band][i] = spectra.R[band][i].dot(mx)

    #- Now combine, if needed, to a single wavelength grid across all cameras
    if spectra.bands == ['brz'] :
        model_wave = spectra.wave['brz']
        mflux = model_flux['brz']

    elif np.all([ band in spectra.bands for band in ['b','r','z'] ]) :
        br_split = 0.5*(spectra.wave['b'][-1] + spectra.wave['r'][0])
        rz_split = 0.5*(spectra.wave['r'][-1] + spectra.wave['z'][0])
        keep = dict()
        keep['b'] = (spectra.wave['b'] < br_split)
        keep['r'] = (br_split <= spectra.wave['r']) & (spectra.wave['r'] < rz_split)
        keep['z'] = (rz_split <= spectra.wave['z'])
        model_wave = np.concatenate( [
            spectra.wave['b'][keep['b']],
            spectra.wave['r'][keep['r']],
            spectra.wave['z'][keep['z']],
        ] )
        mflux = np.concatenate( [
            model_flux['b'][:, keep['b']],
            model_flux['r'][:, keep['r']],
            model_flux['z'][:, keep['z']],
        ], axis=1 )
    else :
        raise RunTimeError("create_model: Set of bands for spectra not supported")

    return model_wave, mflux


def _viewer_urls(spectra, zoom=13, layer='ls-dr9'):
    """Return legacysurvey.org viewer URLs for all spectra.

    Note: `layer` does not apply to the JPEG cutout service.
    """
    u = "https://www.legacysurvey.org/viewer/jpeg-cutout?ra={0:f}&dec={1:f}&zoom={2:d}"
    v = "https://www.legacysurvey.org/viewer/?ra={0:f}&dec={1:f}&zoom={2:d}&layer={3}&mark={0:f},{1:f}"
    if hasattr(spectra, 'fibermap'):
        try:
            ra = spectra.fibermap['RA_TARGET']
            dec = spectra.fibermap['DEC_TARGET']
        except KeyError:
            ra = spectra.fibermap['TARGET_RA']
            dec = spectra.fibermap['TARGET_DEC']
    else:
        ra = spectra.meta['plugmap']['RA']
        dec = spectra.meta['plugmap']['DEC']
    return [(u.format(ra[i], dec[i], zoom, layer),
             v.format(ra[i], dec[i], zoom, layer),
             'RA, Dec = {0:.4f}, {1:+.4f}'.format(ra[i], dec[i]))
            for i in range(len(ra))]


def make_cds_spectra(spectra, with_noise) :
    """ Creates column data source for b,r,z observed spectra """
    cds_spectra = list()
    if isinstance(spectra, SpectrumList):
        s = spectra
        bands = spectra.bands
    else:
        s = [spectra]
        bands = ['coadd']
    for j, band in enumerate(bands):
        cdsdata=dict(origwave=s[j].spectral_axis.value.copy(),
                     plotwave=s[j].spectral_axis.value.copy())
        for i in range(s[j].flux.shape[0]):
            key = 'origflux'+str(i)
            cdsdata[key] = s[j].flux.value[i, :].copy()
            if with_noise :
                key = 'orignoise'+str(i)
                ivar = s[j].uncertainty.array[i, :].copy()
                noise = np.zeros(len(ivar))
                w, = np.where( (ivar > 0))
                noise[w] = 1/np.sqrt(ivar[w])
                cdsdata[key] = noise
        cdsdata['plotflux'] = cdsdata['origflux0']
        if with_noise :
            cdsdata['plotnoise'] = cdsdata['orignoise0']
        cds_spectra.append( bk.ColumnDataSource(cdsdata, name=band) )

    return cds_spectra

def make_cds_coaddcam_spec(spectra, with_noise) :
    """ Creates column data source for camera-coadded observed spectra
        Do NOT store all coadded spectra in CDS obj, to reduce size of html files
        Except for the first spectrum, coaddition is done later in javascript
    """

    coadd_wave, coadd_flux, coadd_ivar = coaddcam_prospect(spectra)
    cds_coaddcam_data = dict(
        origwave = coadd_wave.copy(),
        plotwave = coadd_wave.copy(),
        plotflux = coadd_flux[0,:].copy(),
        plotnoise = np.ones(len(coadd_wave))
    )
    if with_noise :
        w, = np.where( (coadd_ivar[0,:] > 0) )
        cds_coaddcam_data['plotnoise'][w] = 1/np.sqrt(coadd_ivar[0,:][w])
    cds_coaddcam_spec = bk.ColumnDataSource(cds_coaddcam_data)

    return cds_coaddcam_spec

def make_cds_model(model) :
    """ Creates column data source for model spectrum """

    mwave, mflux = model
    cds_model_data = dict(
        origwave = mwave.copy(),
        plotwave = mwave.copy(),
        plotflux = np.zeros(len(mwave)),
    )
    for i in range(len(mflux)):
        key = 'origflux'+str(i)
        cds_model_data[key] = mflux[i]

    cds_model_data['plotflux'] = cds_model_data['origflux0']
    cds_model = bk.ColumnDataSource(cds_model_data)

    return cds_model

def make_cds_targetinfo(spectra, zcatalog, is_coadded, mask_type, username=" ") :
    """ Creates column data source for target-related metadata, from zcatalog, fibermap and VI files """
    target_info = list()
    if isinstance(spectra, Spectrum1D):
        assert mask_type in ['PRIMTARGET', 'SECTARGET',
                             'BOSS_TARGET1', 'BOSS_TARGET2',
                             'ANCILLARY_TARGET1', 'ANCILLARY_TARGET2',
                             'EBOSS_TARGET0', 'EBOSS_TARGET1', 'EBOSS_TARGET2',]
        nspec = spectra.flux.shape[0]
        for i, row in enumerate(spectra.meta['plugmap']):
            target_bit_names = mask_type + ' (DUMMY)'
            target_info.append(target_bit_names)

        cds_targetinfo = bk.ColumnDataSource(
            dict(target_info=target_info),
            name='target_info')

        bands = ['u', 'g', 'r', 'i', 'z']
        for i, bandname in enumerate(bands):
            # mag = np.zeros(len(spectra.meta['plugmap']))
            mag = spectra.meta['plugmap']['MAG'][:, i]
            # extinction = np.ones(len(flux))
            # if ('MW_TRANSMISSION_'+bandname) in spectra.fibermap.keys() :
            #     extinction = spectra.fibermap['MW_TRANSMISSION_'+bandname]
            # w, = np.where( (flux>0) & (extinction>0) )
            # mag[w] = -2.5*np.log10(flux[w]/extinction[w])+22.5
            cds_targetinfo.add(mag, name='mag_'+bandname)

        if zcatalog is not None :
            cds_targetinfo.add(zcatalog['Z'], name='z')
            cds_targetinfo.add(zcatalog['CLASS'].astype('U{0:d}'.format(zcatalog['CLASS'].dtype.itemsize)), name='spectype')
            cds_targetinfo.add(zcatalog['SUBCLASS'].astype('U{0:d}'.format(zcatalog['SUBCLASS'].dtype.itemsize)), name='subtype')
            cds_targetinfo.add(zcatalog['Z_ERR'], name='zerr')
            cds_targetinfo.add(zcatalog['ZWARNING'], name='zwarn')
            cds_targetinfo.add(zcatalog['RCHI2DIFF'], name='deltachi2')
        else :
            cds_targetinfo.add(np.zeros(nspec), name='z')
            cds_targetinfo.add([" " for i in range(nspec)], name='spectype')
            cds_targetinfo.add([" " for i in range(nspec)], name='subtype')
            cds_targetinfo.add(np.zeros(nspec), name='zerr')
            cds_targetinfo.add([0 for i in range(nspec)], name='zwarn')
            cds_targetinfo.add(np.zeros(nspec), name='deltachi2')

        # if not is_coadded and 'EXPID' in spectra.fibermap.keys() :
         #    cds_targetinfo.add(spectra.fibermap['EXPID'], name='expid')
        # else : # If coadd, fill VI accordingly
        cds_targetinfo.add(['-1' for i in range(nspec)], name='expid')
        cds_targetinfo.add([str(x.tolist()) for x in spectra.meta['plugmap']['OBJID']], name='targetid') # !! No int64 in js !!

        #- Get desispec version
        #- TODO : get redrock version (from zcatalog...)
        desispec_specversion = "SDSS"
        # for xx,yy in spectra.meta.items() :
        #     if yy=="desispec" :
        #         desispec_specversion = spectra.meta[xx.replace('NAM','VER')]
        cds_targetinfo.add([desispec_specversion for i in range(nspec)], name='spec_version')
        cds_targetinfo.add(np.zeros(nspec), name='redrock_version')

    else:
        assert mask_type in ['SV1_DESI_TARGET', 'DESI_TARGET', 'CMX_TARGET']
        for i, row in enumerate(spectra.fibermap):
            if mask_type == 'SV1_DESI_TARGET' :
                target_bit_names = ' '.join(sv1_desi_mask.names(row['SV1_DESI_TARGET']))
            elif mask_type == 'DESI_TARGET' :
                target_bit_names = ' '.join(desi_mask.names(row['DESI_TARGET']))
            elif mask_type == 'CMX_TARGET' :
                target_bit_names = ' '.join(cmx_mask.names(row['CMX_TARGET']))
            txt = target_bit_names
            if not is_coadded :
                ## BYPASS DIV
                #           txt += '<BR />'
                if 'NIGHT' in spectra.fibermap.keys() : txt += "Night : {}".format(row['NIGHT'])
                if 'EXPID' in spectra.fibermap.keys() : txt += "Exposure : {}".format(row['EXPID'])
                if 'FIBER' in spectra.fibermap.keys() : txt += "Fiber : {}".format(row['FIBER'])
            target_info.append(txt)

        cds_targetinfo = bk.ColumnDataSource(
            dict(target_info=target_info),
            name='target_info')

        ## BYPASS DIV : Added photometry fields ; also add several bands
        bands = ['G','R','Z', 'W1', 'W2']
        for bandname in bands :
            mag = np.zeros(spectra.num_spectra())
            flux = spectra.fibermap['FLUX_'+bandname]
            extinction = np.ones(len(flux))
            if ('MW_TRANSMISSION_'+bandname) in spectra.fibermap.keys() :
                extinction = spectra.fibermap['MW_TRANSMISSION_'+bandname]
            w, = np.where( (flux>0) & (extinction>0) )
            mag[w] = -2.5*np.log10(flux[w]/extinction[w])+22.5
            cds_targetinfo.add(mag, name='mag_'+bandname)

        nspec = spectra.num_spectra()

        if zcatalog is not None :
            cds_targetinfo.add(zcatalog['Z'], name='z')
            cds_targetinfo.add(zcatalog['SPECTYPE'].astype('U{0:d}'.format(zcatalog['SPECTYPE'].dtype.itemsize)), name='spectype')
            cds_targetinfo.add(zcatalog['SUBTYPE'].astype('U{0:d}'.format(zcatalog['SUBTYPE'].dtype.itemsize)), name='subtype')
            cds_targetinfo.add(zcatalog['ZERR'], name='zerr')
            cds_targetinfo.add(zcatalog['ZWARN'], name='zwarn')
            cds_targetinfo.add(zcatalog['DELTACHI2'], name='deltachi2')
        else :
            cds_targetinfo.add(np.zeros(nspec), name='z')
            cds_targetinfo.add([" " for i in range(nspec)], name='spectype')
            cds_targetinfo.add([" " for i in range(nspec)], name='subtype')
            cds_targetinfo.add(np.zeros(nspec), name='zerr')
            cds_targetinfo.add([0 for i in range(nspec)], name='zwarn')
            cds_targetinfo.add(np.zeros(nspec), name='deltachi2')

        if not is_coadded and 'EXPID' in spectra.fibermap.keys() :
            cds_targetinfo.add(spectra.fibermap['EXPID'], name='expid')
        else : # If coadd, fill VI accordingly
            cds_targetinfo.add(['-1' for i in range(nspec)], name='expid')
        cds_targetinfo.add([str(x) for x in spectra.fibermap['TARGETID']], name='targetid') # !! No int64 in js !!

        #- Get desispec version
        #- TODO : get redrock version (from zcatalog...)
        desispec_specversion = "0"
        for xx,yy in spectra.meta.items() :
            if yy=="desispec" :
                desispec_specversion = spectra.meta[xx.replace('NAM','VER')]
        cds_targetinfo.add([desispec_specversion for i in range(nspec)], name='spec_version')
        cds_targetinfo.add(np.zeros(nspec), name='redrock_version')
        cds_targetinfo.add(np.zeros(nspec)-1, name='template_version')

    # VI inputs
    cds_targetinfo.add([username for i in range(nspec)], name='VI_scanner')
    cds_targetinfo.add(["-1" for i in range(nspec)], name='VI_class_flag')
    cds_targetinfo.add(["" for i in range(nspec)], name='VI_issue_flag')
    cds_targetinfo.add(["" for i in range(nspec)], name='VI_z')
    cds_targetinfo.add(["" for i in range(nspec)], name='VI_spectype')
    cds_targetinfo.add(["" for i in range(nspec)], name='VI_comment')

    return cds_targetinfo


def grid_thumbs(spectra, thumb_width, x_range=(3400,10000), thumb_height=None, resamp_factor=15, ncols_grid=5, titles=None) :
    '''
    Create a bokeh gridplot of thumbnail pictures from spectra
    - coadd arms
    - smooth+resample to reduce size of embedded CDS, according to resamp_factor
    - titles : optional list of titles for each thumb
    '''

    if thumb_height is None : thumb_height = thumb_width//2
    if titles is not None : assert len(titles) == spectra.num_spectra()
    # thumb_wave, thumb_flux, dummy = mycoaddcam(spectra)
    thumb_wave, thumb_flux, dummy = mycoaddcam.coaddcam_prospect(spectra)
    kernel = astropy.convolution.Gaussian1DKernel(stddev=resamp_factor)

    thumb_plots = []
    for i_spec in range(spectra.num_spectra()) :
        x_vals = (thumb_wave[::resamp_factor])[resamp_factor:-resamp_factor]
        # Use astropy convolution : handles NaNs
        y_vals = astropy.convolution.convolve(thumb_flux[i_spec,:], kernel)
        y_vals = (y_vals[::resamp_factor])[resamp_factor:-resamp_factor]
        x_vals = x_vals[~np.isnan(y_vals)] # Needed to avoid 'ValueError: Out of range float values are not JSON compliant'
        y_vals = y_vals[~np.isnan(y_vals)]
        if len(x_vals)==0 : # All NaN ... this should not happen ...
            ymin, ymax = -1, 1
        else :
            yampl = np.max(y_vals) - np.min(y_vals)
            ymin = np.min(y_vals) - 0.1*yampl
            ymax = np.max(y_vals) + 0.1*yampl
        plot_title = None
        if titles is not None : plot_title = titles[i_spec]
        mini_plot = bk.figure(plot_width=thumb_width, plot_height=thumb_height, x_range=x_range, y_range=(ymin,ymax), title=plot_title)
        if len(x_vals)!=0 : mini_plot.line(x_vals, y_vals, line_color='red')
        mini_plot.xaxis.visible = False
        mini_plot.yaxis.visible = False
        mini_plot.min_border_left = 0
        mini_plot.min_border_right = 0
        mini_plot.min_border_top = 0
        mini_plot.min_border_bottom = 0
        thumb_plots.append(mini_plot)

    return bl.gridplot(thumb_plots, ncols=ncols_grid, toolbar_location=None, sizing_mode='scale_width')


def plotspectra(spectra, zcatalog=None, redrock_cat=None, notebook=False, html_dir=None, title=None,
                with_imaging=True, with_noise=True, with_thumb_tab=True,
                with_thumb_only_page=False, with_vi_widgets=True, is_coadded=False,
                with_coaddcam=True, mask_type='DESI_TARGET',
                model_from_zcat=True, model=None, template_dir=None, archetype_fit=False, archetypes_dir=None):
    '''Main prospect routine. From a set of spectra, create a bokeh document
    used for VI, to be displayed as an HTML page or within a Jupyter notebook.

    Parameters
    ----------
    spectra : :class:`~specutils.Spectrum1D` or :class:`~specutils.SpectrumList`
        Input spectra.  :class:`~specutils.SpectrumList` are assumed to be
        DESI spectra.  Otherwise SDSS/BOSS/eBOSS is assumed.
    zcatalog : :class:`~astropy.table.Table`, optional
        Redshift values, matched one-to-one with the input spectra.
    redrock_cat : :class:`~astropy.table.Table`, optional
        Redrock output (as defined in :func:`~prospect.utilities.match_redrock_zfit_to_spectra`).
        Entries must be matched one-by-one (in order) to spectra.
    notebook : :class:`bool`, optional
        If ``True``, bokeh outputs the viewer to a Jupyter notebook.
    html_dir : :class:`str`, optional
        Directory to store the HTML page if `notebook` is ``False``.
    title : :class:`str`, optional
        Title used to name the HTML page / the bokeh figure / the VI file.
    with_imaging : :class:`bool`, optional
        If ``False``, don't include thumb image from http://legacysurvey.org/viewer.
    with_noise : :class:`bool`, optional
        If ``False``, don't include uncertainty for each spectrum.
    with_thumb_tab : :class:`bool`, optional
        If ``False``, don't include a tab with spectra thumbnails.
    with_thumb_only_page : :class:`bool`, optional
        When creating a static HTML (`notebook` is ``False``), a light HTML
        page including only the thumb gallery will also be produced.
    with_vi_widgets : :class:`bool`, optional
        Include widgets used to enter VI information. Set it to ``False`` if
        you do not intend to record VI files.
    is_coadded : :class:`bool`, optional
        Set to ``True`` if `spectra` are coadds.  This will always be assumed
        for SDSS-style inputs, but not for DESI inputs.
    with_coaddcam : :class:`bool`, optional
        Include camera-coaddition, only relevant for DESI.
    mask_type : :class:`str`, optional
        Bitmask type to identify target categories in the spectra. For DESI
        these could be: DESI_TARGET, SV1_DESI_TARGET, CMX_TARGET.
    model_from_zcat : :class:`bool`, optional
        If ``True``, model spectra will be computed from the input `zcatalog`.
    model : :func:`tuple`, optional
        If set, use this input set of model spectra instead of computing it from `zcatalog`.
        model consists of (mwave, mflux); model must be entry-matched to `zcatalog`.
    template_dir : :class:`str`, optional
        Redrock template directory.
    archetype_fit : :class:`bool`, optional
        If ``True``, assume `zcatalog` derived from :command:`redrock --archetypes`
        and plot model accordingly.
    archetypes_dir : :class:`str`, optional
        Directory path for archetypes if not :envvar:`RR_ARCHETYPE_DIR`.
    '''
    # Set masked bins to NaN for compatibility with bokeh.
    if isinstance(spectra, Spectrum1D):
        # We will assume this is from an SDSS/BOSS/eBOSS spPlate file.
        sdss = True
        is_coadded = True
        nspec = spectra.flux.shape[0]
        bad = (spectra.uncertainty.array == 0.0) | spectra.mask
        spectra.flux[bad] = np.nan
    elif isinstance(spectra, SpectrumList):
        # We will assume this is from a DESI spectra-64 file.
        sdss = False
        nspec = spectra[0].flux.shape[0]
        for s in spectra:
            bad = (s.uncertainty.array == 0.0) | s.mask
            s.flux[bad] = np.nan
    else:
        raise ValueError('Unsupported type for input spectra!')

    if title is None:
        title = "specviewer"

    #- Input zcatalog / model
    if zcatalog is not None:
        if sdss:
            if len(zcatalog) != spectra.flux.shape[0]:
                raise ValueError('zcatalog and spectra do not match (different lengths)')
        else:
            if np.any(zcatalog['TARGETID'] != spectra.fibermap['TARGETID']) :
                raise ValueError('zcatalog and spectra do not match (different targetids)')

        if model is not None:
            # SDSS spectra will supply the model.
            assert not model_from_zcat
            mwave, mflux = model
            if len(mflux) != nspec:
                raise ValueError("model fluxes do not match spectra (different nb of entries)")

        if model_from_zcat:
            # DESI spectra will obtain the model from templates.
            model = create_model(spectra, zcatalog,
                                 archetype_fit=archetype_fit,
                                 archetypes_dir=archetypes_dir,
                                 template_dir=template_dir)

    #-----
    #- Initialize Bokeh output
    if notebook:
        assert not with_thumb_only_page
        bk.output_notebook()
    else :
        if html_dir is None : raise RuntimeError("Need html_dir")
        html_page = os.path.join(html_dir, "specviewer_"+title+".html")
        bk.output_file(html_page, title='DESI spectral viewer')

    #-----
    #- Gather information into ColumnDataSource objects for Bokeh
    cds_spectra = make_cds_spectra(spectra, with_noise)
    if with_coaddcam :
        cds_coaddcam_spec = make_cds_coaddcam_spec(spectra, with_noise)
    else :
        cds_coaddcam_spec = None
    if model is not None:
        cds_model = make_cds_model(model)
    else:
        cds_model = None

    if redrock_cat is not None :
        # TODO unhardcode delta_lambd_templates=3
        if np.any(redrock_cat['TARGETID'] != spectra.fibermap['TARGETID']) :
            raise RuntimeError('redrock_cat and spectra do not match (different targetids)')
        if zcatalog is None :
            raise ValueError('Redrock_cat was provided but not zcatalog.')

        with_fit_templates = False if num_approx_fits==0 else True
        template_dicts = make_template_dicts(redrock_cat, delta_lambd_templates=3,
                                             with_fit_templates=with_fit_templates, template_dir=template_dir)
        nfits_redrock_cat = template_dicts[1]['Nfit']
        if num_approx_fits is None : num_approx_fits = nfits_redrock_cat
        if (num_approx_fits > nfits_redrock_cat) : raise ValueError("num_approx_fits too large wrt redrock_cat")
        if with_full_2ndfit :
            zcat_2ndfit = create_zcat_from_redrock_cat(redrock_cat, fit_num=1)
            model_2ndfit = create_model(spectra, zcat_2ndfit, archetype_fit=archetype_fit,
                                        archetypes_dir=archetypes_dir, template_dir=template_dir)
            cds_model_2ndfit = make_cds_model(model_2ndfit)
        else :
            cds_model_2ndfit = None
        # Now the "plot" CDS : initialize it to best fit
        cds_othermodel = bk.ColumnDataSource({
            'plotwave' : cds_model.data['plotwave'],
            'origwave' : cds_model.data['origwave'],
            'origflux' : cds_model.data['origflux0'],
            'plotflux' : cds_model.data['origflux0'],
            'zref' : zcatalog['Z'][0]+np.zeros(len(cds_model.data['origflux0'])) # trick to track the z reference in model
        })
    else :
        cds_model_2ndfit = None
        template_dicts = None
        cds_othermodel =  None

    if notebook and ("USER" in os.environ) :
        username = os.environ['USER'][0:3] # 3-letter acronym
    else :
        username = " "
    cds_targetinfo = make_cds_targetinfo(spectra, zcatalog, is_coadded, mask_type, username=username)

    #-------------------------
    #-- Graphical objects --
    #-------------------------


    #-----
    #- Main figure
    #- Determine initial ymin, ymax, xmin, xmax
    ymin = ymax = xmax = 0.0
    xmin = 100000.
    xmargin = 300.
    if sdss:
        bands = ['coadd']
        ymin = np.nanmin(spectra.flux.value[0])
        ymax = np.nanmax(spectra.flux.value[0])
        xmin = np.min(spectra.spectral_axis.value)
        xmax = np.max(spectra.spectral_axis.value)
    else:
        bands = spectra.bands
        for i, band in enumerate(bands):
            ymin = min(ymin, np.nanmin(spectra[i].flux.value[0]))
            ymax = max(ymax, np.nanmax(spectra[i].flux.value[0]))
            xmin = min(xmin, np.min(spectra[i].spectral_axis.value))
            xmax = max(xmax, np.max(spectra[i].spectral_axis.value))
    xmin -= xmargin
    xmax += xmargin

    plot_width=800
    plot_height=400
    plot_widget_width = (plot_width+(plot_height//2))//2 - 40 # used for widgets scaling
    tools = 'pan,box_zoom,wheel_zoom,save'
    tooltips_fig = [("wave","$x"),("flux","$y")]
    fig = bk.figure(height=plot_height, width=plot_width, title=title,
        tools=tools, toolbar_location='above', tooltips=tooltips_fig, y_range=(ymin, ymax), x_range=(xmin, xmax))
    fig.sizing_mode = 'stretch_width'
    fig.toolbar.active_drag = fig.tools[0]    #- pan zoom (previously box)
    fig.toolbar.active_scroll = fig.tools[2]  #- wheel zoom
    fig.xaxis.axis_label = 'Wavelength [Å]'
    fig.yaxis.axis_label = 'Flux'
    fig.xaxis.axis_label_text_font_style = 'normal'
    fig.yaxis.axis_label_text_font_style = 'normal'
    colors = dict(b='#1f77b4', r='#d62728', z='maroon', coadd='#d62728', brz='#d62728')
    noise_colors = dict(b='greenyellow', r='green', z='forestgreen', coadd='green', brz='green')
    alpha_discrete = 0.2 # alpha for "almost-hidden" curves (single-arm spectra and noise by default)
    if not with_coaddcam : alpha_discrete = 1

    #- Highlight overlap regions between arms
    ## overlap wavelengths are hardcoded, from 1907.10688 (Table 1)
    overlap_waves = [ [5660, 5930], [7470, 7720] ]
    alpha_overlapband = 0.03
    overlap_bands = []
    if bands == ['brz'] :
        for i in range(len(overlap_waves)) :
            fill_alpha = alpha_overlapband if with_coaddcam else 0
            overlap_bands.append( BoxAnnotation(left=overlap_waves[i][0], right=overlap_waves[i][1], fill_color='blue', fill_alpha=fill_alpha, line_alpha=0) )
            fig.add_layout(overlap_bands[-1])

    data_lines = list()
    for spec in cds_spectra:
        lx = fig.line('plotwave', 'plotflux', source=spec, line_color=colors[spec.name], line_alpha=alpha_discrete)
        data_lines.append(lx)
    if with_coaddcam :
        lx = fig.line('plotwave', 'plotflux', source=cds_coaddcam_spec, line_color=colors['coadd'], line_alpha=1)
        data_lines.append(lx)

    noise_lines = list()
    if with_noise :
        for spec in cds_spectra :
            lx = fig.line('plotwave', 'plotnoise', source=spec, line_color=noise_colors[spec.name], line_alpha=alpha_discrete)
            noise_lines.append(lx)
        if with_coaddcam :
            lx = fig.line('plotwave', 'plotnoise', source=cds_coaddcam_spec, line_color=noise_colors['coadd'], line_alpha=1)
            noise_lines.append(lx)

    model_lines = list()
    if cds_model is not None:
        lx = fig.line('plotwave', 'plotflux', source=cds_model, line_color='black')
        model_lines.append(lx)

    othermodel_lines = list()
    if cds_othermodel is not None :
        lx = fig.line('plotwave', 'plotflux', source=cds_othermodel, line_color='black', line_dash='dashed')
        othermodel_lines.append(lx)

    legend_items = [("data",  data_lines[-1::-1])] #- reversed to get blue as lengend entry
    if cds_model is not None :
        legend_items.append(("pipeline fit", model_lines))
    if cds_othermodel is not None :
        legend_items.append(("other model", othermodel_lines))
    if with_noise :
        legend_items.append(("noise", noise_lines[-1::-1])) # same as for data_lines
    legend = Legend(items=legend_items)
    fig.add_layout(legend, 'center')
    fig.legend.click_policy = 'hide'    #- or 'mute'

    #-----
    #- Zoom figure around mouse hover of main plot
    tooltips_zoomfig = [("wave","$x"),("flux","$y")]
    zoomfig = bk.figure(height=plot_height//2, width=plot_height//2,
        y_range=fig.y_range, x_range=(5000,5100),
        # output_backend="webgl",
        toolbar_location=None, tooltips=tooltips_zoomfig, tools=[])

    zoom_data_lines = list()
    zoom_noise_lines = list()
    for spec in cds_spectra:
        zoom_data_lines.append(zoomfig.line('plotwave', 'plotflux', source=spec,
            line_color=colors[spec.name], line_width=1, line_alpha=alpha_discrete))
        if with_noise :
            zoom_noise_lines.append(zoomfig.line('plotwave', 'plotnoise', source=spec,
                            line_color=noise_colors[spec.name], line_width=1, line_alpha=alpha_discrete))
    if with_coaddcam :
        zoom_data_lines.append(zoomfig.line('plotwave', 'plotflux', source=cds_coaddcam_spec, line_color=colors['coadd'], line_alpha=1))
        if with_noise :
            lx = zoomfig.line('plotwave', 'plotnoise', source=cds_coaddcam_spec, line_color=noise_colors['coadd'], line_alpha=1)
            zoom_noise_lines.append(lx)

    if cds_model is not None:
        lx = zoomfig.line('plotwave', 'plotflux', source=cds_model, line_color='black')
    if cds_othermodel is not None :
        lx = zoomfig.line('plotwave', 'plotflux', source=cds_othermodel, line_color='black', line_dash='dashed')

    #- Callback to update zoom window x-range
    zoom_callback = CustomJS(
        args=dict(zoomfig=zoomfig,fig=fig),
        code="""
            zoomfig.x_range.start = cb_obj.x - 100;
            zoomfig.x_range.end = cb_obj.x + 100;
        """)

    fig.js_on_event(bokeh.events.MouseMove, zoom_callback)

    #-----
    #- Targeting image
    if with_imaging :
        imfig = bk.figure(width=plot_height//2, height=plot_height//2,
                          x_range=(0, 256), y_range=(0, 256),
                          x_axis_location=None, y_axis_location=None,
                          output_backend="webgl",
                          toolbar_location=None, tools=[])
        imfig.min_border_left = 0
        imfig.min_border_right = 0
        imfig.min_border_top = 0
        imfig.min_border_bottom = 0

        imfig_urls = _viewer_urls(spectra)
        imfig_source = ColumnDataSource(data=dict(url=[imfig_urls[0][0]],
                                                  txt=[imfig_urls[0][2]]))

        imfig_img = imfig.image_url('url', source=imfig_source, x=1, y=1, w=256, h=256, anchor='bottom_left')
        imfig_txt = imfig.text(10, 256-30, text='txt', source=imfig_source,
                               text_color='yellow', text_font_size='8pt')
        # cross-hair
        imfig.multi_line([[129-15,129-5],[129+15,129+5],[129,129],[129,129]],
                         [[129,129],[129,129],[129-15,129-5],[129+5,129+15]], line_width=1, line_color='yellow')
    else :
        imfig = bl.Spacer(width=plot_height//2, height=plot_height//2)
        imfig_source = imfig_urls = None

    #-----
    #- Emission and absorption lines
    z = zcatalog['Z'][0] if (zcatalog is not None) else 0.0
    line_data, lines, line_labels = add_lines(fig, z=z)
    zoom_line_data, zoom_lines, zoom_line_labels = add_lines(zoomfig, z=z, label_offsets=[50, 5])


    #-------------------------
    #-- Widgets and callbacks --
    #-------------------------

    js_files = get_resources('js')

    #-----
    #- Ifiberslider and smoothing widgets
    # Ifiberslider's value controls which spectrum is displayed
    # These two widgets call update_plot(), later defined
    slider_end = nspec-1 if nspec > 1 else 0.5 # Slider cannot have start=end
    ifiberslider = Slider(start=0, end=slider_end, value=0, step=1, title='Spectrum (of '+str(nspec)+')')
    smootherslider = Slider(start=0, end=26, value=0, step=1.0, title='Gaussian Sigma Smooth')

    #-----
    #- Navigation buttons
    navigation_button_width = 30
    prev_button = Button(label="<", width=navigation_button_width)
    next_button = Button(label=">", width=navigation_button_width)
    prev_callback = CustomJS(
        args=dict(ifiberslider=ifiberslider),
        code="""
        if(ifiberslider.value>0 && ifiberslider.end>=1) {
            ifiberslider.value--
        }
        """)
    next_callback = CustomJS(
        args=dict(ifiberslider=ifiberslider, nspec=nspec),
        code="""
        if(ifiberslider.value<nspec-1 && ifiberslider.end>=1) {
            ifiberslider.value++
        }
        """)
    prev_button.js_on_event('button_click', prev_callback)
    next_button.js_on_event('button_click', next_callback)


    #-----
    #- Axis reset button (superseeds the default bokeh "reset"
    reset_plotrange_button = Button(label="Reset X-Y range", button_type="default")
    reset_plotrange_code = js_files["adapt_plotrange.js"] + js_files["reset_plotrange.js"]
    reset_plotrange_callback = CustomJS(args = dict(fig=fig, xmin=xmin, xmax=xmax, spectra=cds_spectra),
                                        code = reset_plotrange_code)
    reset_plotrange_button.js_on_event('button_click', reset_plotrange_callback)


    #-----
    #- Redshift / wavelength scale widgets
    z1 = np.floor(z*100)/100
    dz = z-z1
    zslider = Slider(start=-0.1, end=5.0, value=z1, step=0.01, title='Redshift rough tuning')
    dzslider = Slider(start=0.0, end=0.0099, value=dz, step=0.0001, title='Redshift fine-tuning')
    dzslider.format = "0[.]0000"
    z_input = TextInput(value="{:.4f}".format(z), title="Redshift value:")

    #- Observer vs. Rest frame wavelengths
    waveframe_buttons = RadioButtonGroup(
        labels=["Obs", "Rest"], active=0)

    zslider_callback  = CustomJS(
        args=dict(zslider=zslider, dzslider=dzslider, z_input=z_input),
        code="""
        // Protect against 1) recursive call with z_input callback;
        //   2) out-of-range zslider values (should never happen in principle)
        var z1 = Math.floor(parseFloat(z_input.value)*100) / 100
        if ( (Math.abs(zslider.value-z1) >= 0.01) &&
             (zslider.value >= -0.1) && (zslider.value <= 5.0) ){
             var new_z = zslider.value + dzslider.value
             z_input.value = new_z.toFixed(4)
            }
        """)

    dzslider_callback  = CustomJS(
        args=dict(zslider=zslider, dzslider=dzslider, z_input=z_input),
        code="""
        var z = parseFloat(z_input.value)
        var z1 = Math.floor(z) / 100
        var z2 = z-z1
        if ( (Math.abs(dzslider.value-z2) >= 0.0001) &&
             (dzslider.value >= 0.0) && (dzslider.value <= 0.0099) ){
             var new_z = zslider.value + dzslider.value
             z_input.value = new_z.toFixed(4)
            }
        """)

    zslider.js_on_change('value', zslider_callback)
    dzslider.js_on_change('value', dzslider_callback)

    z_button_width = 30
    z_minus_button = Button(label="<", width=z_button_width)
    z_plus_button = Button(label=">", width=z_button_width)
    z_minus_callback = CustomJS(
        args=dict(z_input=z_input),
        code="""
        var z = parseFloat(z_input.value)
        if(z >= -0.09) {
            z -= 0.01
            z_input.value = z.toFixed(4)
        }
        """)
    z_plus_callback = CustomJS(
        args=dict(z_input=z_input),
        code="""
        var z = parseFloat(z_input.value)
        if(z <= 4.99) {
            z += 0.01
            z_input.value = z.toFixed(4)
        }
        """)
    z_minus_button.js_on_event('button_click', z_minus_callback)
    z_plus_button.js_on_event('button_click', z_plus_callback)

    zreset_button = Button(label='Reset to z_pipe')
    zreset_callback = CustomJS(
        args=dict(z_input=z_input, targetinfo=cds_targetinfo, ifiberslider=ifiberslider),
        code="""
            var ifiber = ifiberslider.value
            var z = targetinfo.data['z'][ifiber]
            z_input.value = z.toFixed(4)
        """)
    zreset_button.js_on_event('button_click', zreset_callback)

    z_input_callback = CustomJS(
        args=dict(spectra = cds_spectra,
            coaddcam_spec = cds_coaddcam_spec,
            model = cds_model,
            othermodel = cds_othermodel,
            targetinfo = cds_targetinfo,
            ifiberslider = ifiberslider,
            zslider=zslider,
            dzslider=dzslider,
            z_input=z_input,
            waveframe_buttons=waveframe_buttons,
            line_data=line_data, lines=lines, line_labels=line_labels,
            zlines=zoom_lines, zline_labels=zoom_line_labels,
            overlap_waves=overlap_waves, overlap_bands=overlap_bands,
            fig=fig
            ),
        code="""
            var z = parseFloat(z_input.value)
            if ( z >=-0.1 && z <= 5.0 ) {
                // update zsliders only if needed (avoid recursive call)
                z_input.value = parseFloat(z_input.value).toFixed(4)
                var z1 = Math.floor(z*100) / 100
                var z2 = z-z1
                if ( Math.abs(z1-zslider.value) >= 0.01) zslider.value = parseFloat(parseFloat(z1).toFixed(2))
                if ( Math.abs(z2-dzslider.value) >= 0.0001) dzslider.value = parseFloat(parseFloat(z2).toFixed(4))
            } else {
                if (z_input.value < -0.1) z_input.value = (-0.1).toFixed(4)
                if (z_input.value > 5) z_input.value = (5.0).toFixed(4)
            }

            var line_restwave = line_data.data['restwave']
            var ifiber = ifiberslider.value
            var waveshift_lines = (waveframe_buttons.active == 0) ? 1+z : 1 ;
            var waveshift_spec = (waveframe_buttons.active == 0) ? 1 : 1/(1+z) ;

            for(var i=0; i<line_restwave.length; i++) {
                lines[i].location = line_restwave[i] * waveshift_lines
                line_labels[i].x = line_restwave[i] * waveshift_lines
                zlines[i].location = line_restwave[i] * waveshift_lines
                zline_labels[i].x = line_restwave[i] * waveshift_lines
            }
            if (overlap_bands.length>0) {
                for (var i=0; i<overlap_bands.length; i++) {
                    overlap_bands[i].left = overlap_waves[i][0] * waveshift_spec
                    overlap_bands[i].right = overlap_waves[i][1] * waveshift_spec
                }
            }

            function shift_plotwave(cds_spec, waveshift) {
                var data = cds_spec.data
                var origwave = data['origwave']
                var plotwave = data['plotwave']
                if ( plotwave[0] != origwave[0] * waveshift ) { // Avoid redo calculation if not needed
                    for (var j=0; j<plotwave.length; j++) {
                        plotwave[j] = origwave[j] * waveshift ;
                    }
                    cds_spec.change.emit()
                }
            }

            for(var i=0; i<spectra.length; i++) {
                shift_plotwave(spectra[i], waveshift_spec)
            }
            if (coaddcam_spec) shift_plotwave(coaddcam_spec, waveshift_spec)

            // Update model wavelength array
            // NEW : don't shift model if othermodel is there
            if (othermodel) {
                var zref = othermodel.data['zref'][0]
                var waveshift_model = (waveframe_buttons.active == 0) ? (1+z)/(1+zref) : 1/(1+zref) ;
                shift_plotwave(othermodel, waveshift_model)
            } else if (model) {
                var zfit = 0.0
                if(targetinfo.data['z'] != undefined) {
                    zfit = targetinfo.data['z'][ifiber]
                }
                var waveshift_model = (waveframe_buttons.active == 0) ? (1+z)/(1+zfit) : 1/(1+zfit) ;
                shift_plotwave(model, waveshift_model)
            }
        """)
    z_input.js_on_change('value', z_input_callback)
    waveframe_buttons.js_on_click(z_input_callback)

    plotrange_callback = CustomJS(
        args = dict(
            z_input=z_input,
            waveframe_buttons=waveframe_buttons,
            fig=fig,
        ),
        code="""
        var z =  parseFloat(z_input.value)
        // Observer Frame
        if(waveframe_buttons.active == 0) {
            fig.x_range.start = fig.x_range.start * (1+z)
            fig.x_range.end = fig.x_range.end * (1+z)
        } else {
            fig.x_range.start = fig.x_range.start / (1+z)
            fig.x_range.end = fig.x_range.end / (1+z)
        }
        """
    )
    waveframe_buttons.js_on_click(plotrange_callback)


    #------
    #- Zoom on the OII doublet TODO mv js code to other file
    # TODO: is there another trick than using a cds to pass the "oii_saveinfo" ?
    # TODO: optimize smoothing for autozoom (current value: 0)
    cds_oii_saveinfo = bk.ColumnDataSource(
        {'xmin':[fig.x_range.start], 'xmax':[fig.x_range.end], 'nsmooth':[smootherslider.value]})
    oii_zoom_button = Button(label="OII-zoom", button_type="default")
    oii_zoom_callback = CustomJS(
        args = dict(z_input=z_input, fig=fig, smootherslider=smootherslider,
                   cds_oii_saveinfo=cds_oii_saveinfo),
        code = """
        // Save previous setting (for the "Undo" button)
        cds_oii_saveinfo.data['xmin'] = [fig.x_range.start]
        cds_oii_saveinfo.data['xmax'] = [fig.x_range.end]
        cds_oii_saveinfo.data['nsmooth'] = [smootherslider.value]
        // Center on the middle of the redshifted OII doublet (vaccum)
        var z = parseFloat(z_input.value)
        fig.x_range.start = 3728.48 * (1+z) - 30
        fig.x_range.end = 3728.48 * (1+z) + 30
        // No smoothing (this implies a call to update_plot)
        smootherslider.value = 0
        """)
    oii_zoom_button.js_on_event('button_click', oii_zoom_callback)

    oii_undo_button = Button(label="Undo", button_type="default")
    oii_undo_callback = CustomJS(
        args = dict(fig=fig, smootherslider=smootherslider, cds_oii_saveinfo=cds_oii_saveinfo),
        code = """
        fig.x_range.start = cds_oii_saveinfo.data['xmin'][0]
        fig.x_range.end = cds_oii_saveinfo.data['xmax'][0]
        smootherslider.value = cds_oii_saveinfo.data['nsmooth'][0]
        """)
    oii_undo_button.js_on_event('button_click', oii_undo_callback)

    #-----
    #- Targeting image callback
    if with_imaging :
        imfig_callback = CustomJS(args=dict(urls=imfig_urls,
                                            ifiberslider=ifiberslider),
                                  code='''window.open(urls[ifiberslider.value][1], "_blank");''')
        imfig.js_on_event('tap', imfig_callback)


    #-----
    #- Highlight individual-arm or camera-coadded spectra
    coaddcam_labels = []
    if cds_coaddcam_spec is not None : coaddcam_labels = ["Camera-coadded", "Single-arm"]
    coaddcam_buttons = RadioButtonGroup(labels=coaddcam_labels, active=0)
    coaddcam_callback = CustomJS(
        args = dict(coaddcam_buttons=coaddcam_buttons,
                    list_lines=[data_lines, noise_lines, zoom_data_lines, zoom_noise_lines],
                    alpha_discrete=alpha_discrete,
                    overlap_bands=overlap_bands,
                    alpha_overlapband=alpha_overlapband),
        code="""
        var n_lines = list_lines[0].length
        for (var i=0; i<n_lines; i++) {
            var new_alpha = 1
            if (coaddcam_buttons.active == 0 && i<n_lines-1) new_alpha = alpha_discrete
            if (coaddcam_buttons.active == 1 && i==n_lines-1) new_alpha = alpha_discrete
            for (var j=0; j<list_lines.length; j++) {
                list_lines[j][i].glyph.line_alpha = new_alpha
            }
        }
        var new_alpha = 0
        if (coaddcam_buttons.active == 0) new_alpha = alpha_overlapband
        for (var j=0; j<overlap_bands.length; j++) {
                overlap_bands[j].fill_alpha = new_alpha
        }
        """
    )
    coaddcam_buttons.js_on_click(coaddcam_callback)

    #-----
    # Display object-related informations
    ## BYPASS DIV to be able to copy targetid...
    ## target_info_div = Div(text=cds_targetinfo.data['target_info'][0])
    tmp_dict = dict()
    tmp_dict['TARGETID'] = [ cds_targetinfo.data['targetid'][0] ]
    tmp_dict['Target class'] = [ cds_targetinfo.data['target_info'][0] ]
    targ_disp_cols = [ TableColumn(field='TARGETID', title='TARGETID', width=150),
                     TableColumn(field='Target class', title='Target class', width=250) ] # TODO tune width
    if sdss:
        phot_bands = ['u', 'g', 'r', 'i', 'z']
    else:
        phot_bands = ['G', 'R', 'Z', 'W1', 'W2']
    for band in phot_bands:
        tmp_dict['mag_'+band] = [ "{:.2f}".format(cds_targetinfo.data['mag_'+band][0]) ]
        targ_disp_cols.append( TableColumn(field='mag_'+band, title='mag_'+band, width=40) )
    targ_disp_cds = bk.ColumnDataSource(tmp_dict, name='targ_disp_cds')
    targ_display = DataTable(source = targ_disp_cds, columns=targ_disp_cols,index_position=None, selectable=True, editable=True) # width=...
    targ_display.height = 2 * targ_display.row_height
    if zcatalog is not None :
        if template_dicts is not None : # Add other best fits
            fit_results = template_dicts[1]
            # Case of DeltaChi2 : compute it from Chi2s
            #    The "DeltaChi2" in rr fits is between best fits for a given (spectype,subtype)
            #    Convention: DeltaChi2 = -1 for the last fit.
            chi2s = fit_results['CHI2'][0]
            full_deltachi2s = np.zeros(len(chi2s))-1
            full_deltachi2s[:-1] = chi2s[1:]-chi2s[:-1]
            tmp_dict = dict(Nfit = np.arange(1,len(chi2s)+1),
                            SPECTYPE = fit_results['SPECTYPE'][0],  # [0:num_best_fits] (if we want to restrict... TODO?)
                            SUBTYPE = fit_results['SUBTYPE'][0],
                            Z = [ "{:.4f}".format(x) for x in fit_results['Z'][0] ],
                            ZERR = [ "{:.4f}".format(x) for x in fit_results['ZERR'][0] ],
                            ZWARN = fit_results['ZWARN'][0],
                            CHI2 = [ "{:.1f}".format(x) for x in fit_results['CHI2'][0] ],
                            DeltaChi2 = [ "{:.1f}".format(x) for x in full_deltachi2s ])
        else :
            tmp_dict = dict(SPECTYPE = [ cds_targetinfo.data['spectype'][0] ],
                SUBTYPE = [ cds_targetinfo.data['subtype'][0] ],
                Z = [ "{:.4f}".format(cds_targetinfo.data['z'][0]) ],
                ZERR = [ "{:.4f}".format(cds_targetinfo.data['zerr'][0]) ],
                ZWARN = [ cds_targetinfo.data['zwarn'][0] ],
                DeltaChi2 = [ "{:.1f}".format(cds_targetinfo.data['deltachi2'][0]) ])
        zcat_disp_cds = bk.ColumnDataSource(tmp_dict, name='zcat_disp_cds')
        zcat_disp_cols = [ TableColumn(field=x, title=t, width=w) for x,t,w in [ ('SPECTYPE','SPECTYPE',70), ('SUBTYPE','SUBTYPE',60), ('Z','Z',50) , ('ZERR','ZERR',50), ('ZWARN','ZWARN',50), ('DeltaChi2','Δχ2(N/N+1)',70)] ]
        if template_dicts is not None :
            zcat_disp_cols.insert(0, TableColumn(field='Nfit', title='Nfit', width=5))
        zcat_display = DataTable(source=zcat_disp_cds, columns=zcat_disp_cols, selectable=False, index_position=None, width=plot_widget_width)
        zcat_display.height = 2 * zcat_display.row_height
        if template_dicts is not None : zcat_display.height = 3 * zcat_display.row_height
    else :
        zcat_display = Div(text="Not available ")
        zcat_disp_cds = None

    #-----
    #- Toggle lines
    lines_button_group = CheckboxButtonGroup(
            labels=["Emission lines", "Absorption lines"], active=[])
    majorline_checkbox = CheckboxGroup(
            labels=['Show only major lines'], active=[])

    lines_callback = CustomJS(
        args = dict(line_data=line_data, lines=lines, line_labels=line_labels, zlines=zoom_lines,
                    zline_labels=zoom_line_labels, lines_button_group=lines_button_group, majorline_checkbox=majorline_checkbox),
        code="""
        var show_emission = false
        var show_absorption = false
        if (lines_button_group.active.indexOf(0) >= 0) {  // index 0=Emission in active list
            show_emission = true
        }
        if (lines_button_group.active.indexOf(1) >= 0) {  // index 1=Absorption in active list
            show_absorption = true
        }

        for(var i=0; i<lines.length; i++) {
            if ( !(line_data.data['major'][i]) && (majorline_checkbox.active.indexOf(0)>=0) ) {
                lines[i].visible = false
                line_labels[i].visible = false
                zlines[i].visible = false
                zline_labels[i].visible = false
            } else if (line_data.data['emission'][i]) {
                lines[i].visible = show_emission
                line_labels[i].visible = show_emission
                zlines[i].visible = show_emission
                zline_labels[i].visible = show_emission
            } else {
                lines[i].visible = show_absorption
                line_labels[i].visible = show_absorption
                zlines[i].visible = show_absorption
                zline_labels[i].visible = show_absorption
            }
        }
        """
    )
    lines_button_group.js_on_click(lines_callback)
    majorline_checkbox.js_on_click(lines_callback)

    #------
    #- Select secondary model to display
    if template_dicts is not None :
        model_options = ['Best fit', '2nd best fit']
        for i in range(1,1+num_approx_fits) :
            ith = 'th'
            if i==1 : ith='st'
            if i==2 : ith='nd'
            if i==3 : ith='rd'
            model_options.append(str(i)+ith+' fit (approx)')
        if with_full_2ndfit is False :
            model_options.remove('2nd best fit')
        for std_template in ['QSO', 'GALAXY', 'STAR'] :
            model_options.append('STD '+std_template)
        model_select = Select(value=model_options[0], title="Other model (dashed curve):", options=model_options)
        cds_median_spectra = make_cds_median_spectra(spectra)
        model_select_code = js_files["interp_grid.js"] + js_files["smooth_data.js"] + js_files["select_model.js"]
        model_select_callback = CustomJS(
            args=dict(ifiberslider = ifiberslider,
                      model_select = model_select,
                      fit_templates=template_dicts[0],
                      cds_othermodel=cds_othermodel,
                      cds_model_2ndfit=cds_model_2ndfit,
                      cds_model = cds_model,
                      fit_results=template_dicts[1],
                      std_templates=template_dicts[2],
                      median_spectra = cds_median_spectra,
                      smootherslider = smootherslider,
                      z_input = z_input,
                      cds_targetinfo = cds_targetinfo),
            code=model_select_code)
        model_select.js_on_change('value',model_select_callback)
    else :
        model_select = None


    #-----
    #- VI-related widgets

    vi_class_labels = [ x["label"] for x in vi_flags if x["type"]=="class" ]
    vi_issue_labels = [ x["label"] for x in vi_flags if x["type"]=="issue" ]
    vi_issue_slabels = [ x["shortlabel"] for x in vi_flags if x["type"]=="issue" ]

    #- VI file name
    default_vi_filename = "desi-vi_"+title
    if username.strip()!="" :
        default_vi_filename += ("_"+username)
    else :
        default_vi_filename += "_unknown-user"
    default_vi_filename += ".csv"
    vi_filename_input = TextInput(value=default_vi_filename, title="VI file name:")

    #- Optional VI flags (issues)
    vi_issue_input = CheckboxGroup(labels=vi_issue_labels, active=[])
    vi_issue_code = js_files["CSVtoArray.js"] + js_files["save_vi.js"]
    vi_issue_code += """
        var issues = []
        for (var i=0; i<vi_issue_labels.length; i++) {
            if (vi_issue_input.active.indexOf(i) >= 0) issues.push(vi_issue_slabels[i])
        }
        if (issues.length > 0) {
            cds_targetinfo.data['VI_issue_flag'][ifiberslider.value] = ( issues.join('') )
        } else {
            cds_targetinfo.data['VI_issue_flag'][ifiberslider.value] = " "
        }
        autosave_vi_localStorage(vi_file_fields, cds_targetinfo.data, title)
        cds_targetinfo.change.emit()
        """
    vi_issue_callback = CustomJS(
        args=dict(cds_targetinfo=cds_targetinfo,ifiberslider = ifiberslider,
                vi_issue_input=vi_issue_input, vi_issue_labels=vi_issue_labels,
                vi_issue_slabels=vi_issue_slabels,
                title=title, vi_file_fields = vi_file_fields),
        code=vi_issue_code )
    vi_issue_input.js_on_click(vi_issue_callback)

    #- Optional VI information on redshift
    vi_z_input = TextInput(value='', title="VI redshift:")
    vi_z_code = js_files["CSVtoArray.js"] + js_files["save_vi.js"]
    vi_z_code += """
        cds_targetinfo.data['VI_z'][ifiberslider.value]=vi_z_input.value
        autosave_vi_localStorage(vi_file_fields, cds_targetinfo.data, title)
        cds_targetinfo.change.emit()
        """
    vi_z_callback = CustomJS(
        args=dict(cds_targetinfo=cds_targetinfo, ifiberslider = ifiberslider, vi_z_input=vi_z_input,
                  title=title, vi_file_fields=vi_file_fields),
        code=vi_z_code )
    vi_z_input.js_on_change('value',vi_z_callback)

    # Copy z value from redshift slider to VI
    z_tovi_button = Button(label='Copy z to VI')
    z_tovi_callback = CustomJS(
        args=dict(z_input=z_input, vi_z_input=vi_z_input),
        code="""
            vi_z_input.value = z_input.value
        """)
    z_tovi_button.js_on_event('button_click', z_tovi_callback)

    #- Optional VI information on spectral type
    vi_category_select = Select(value=' ', title="VI spectype:", options=([' '] + vi_spectypes))
    # The default value set to ' ' as setting value='' does not seem to work well with Select.
    vi_category_code = js_files["CSVtoArray.js"] + js_files["save_vi.js"]
    vi_category_code += """
        if (vi_category_select.value == ' ') {
            cds_targetinfo.data['VI_spectype'][ifiberslider.value]=''
        } else {
            cds_targetinfo.data['VI_spectype'][ifiberslider.value]=vi_category_select.value
        }
        autosave_vi_localStorage(vi_file_fields, cds_targetinfo.data, title)
        cds_targetinfo.change.emit()
        """
    vi_category_callback = CustomJS(
        args=dict(cds_targetinfo=cds_targetinfo, ifiberslider = ifiberslider,
                  vi_category_select=vi_category_select,
                  title=title, vi_file_fields=vi_file_fields),
        code=vi_category_code )
    vi_category_select.js_on_change('value',vi_category_callback)

    #- Optional VI comment
    vi_comment_input = TextInput(value='', title="VI comment (see guidelines):")
    vi_comment_code = js_files["CSVtoArray.js"] + js_files["save_vi.js"]
    vi_comment_code += """
        var stored_comment = (vi_comment_input.value).replace(/./g, function(char){
            if ( char==',' ) {
                return ';'
            } else if ( char.charCodeAt(0)<=127 ) {
                return char
            } else {
                var char_list = ['Å','α','β','γ','δ','λ']
                var char_replace = ['Angstrom','alpha','beta','gamma','delta','lambda']
                for (var i=0; i<char_list.length; i++) {
                    if ( char==char_list[i] ) return char_replace[i]
                }
                return '?'
            }
        })
        cds_targetinfo.data['VI_comment'][ifiberslider.value] = stored_comment
        autosave_vi_localStorage(vi_file_fields, cds_targetinfo.data, title)
        cds_targetinfo.change.emit()
        """
    vi_comment_callback = CustomJS(
        args=dict(cds_targetinfo=cds_targetinfo, ifiberslider = ifiberslider, vi_comment_input=vi_comment_input,
                  title=title, vi_file_fields=vi_file_fields),
        code=vi_comment_code )
    vi_comment_input.js_on_change('value',vi_comment_callback)

    #- List of "standard" VI comment
    vi_std_comment_select = Select(value=" ", title="Standard comment:", options=([' '] + vi_std_comments))
    vi_std_comment_code = """
        if (vi_std_comment_select.value != ' ') {
            if (vi_comment_input.value != '') {
                vi_comment_input.value = vi_comment_input.value + " " + vi_std_comment_select.value
            } else {
                vi_comment_input.value = vi_std_comment_select.value
            }
        }
        """
    vi_std_comment_callback = CustomJS(
        args = dict(vi_std_comment_select=vi_std_comment_select, vi_comment_input=vi_comment_input),
        code = vi_std_comment_code )
    vi_std_comment_select.js_on_change('value', vi_std_comment_callback)

    #- Main VI classification
    vi_class_input = RadioButtonGroup(labels=vi_class_labels)
    vi_class_code = js_files["CSVtoArray.js"] + js_files["save_vi.js"]
    vi_class_code += """
        if ( vi_class_input.active >= 0 ) {
            cds_targetinfo.data['VI_class_flag'][ifiberslider.value] = vi_class_labels[vi_class_input.active]
            //if ( vi_class_labels[vi_class_input.active]=="4" && model) { // Flag '4' => VI_z = z_pipe (if available)
            //    var z = targetinfo.data['z'][ifiberslider.value]
            //    vi_z_input.value = parseFloat(z).toFixed(4)
            //    vi_category_select.value = targetinfo.data['spectype'][ifiberslider.value]
            //}
        } else {
            cds_targetinfo.data['VI_class_flag'][ifiberslider.value] = "-1"
        }
        autosave_vi_localStorage(vi_file_fields, cds_targetinfo.data, title)
        cds_targetinfo.change.emit()
    """
    vi_class_callback = CustomJS(
        args=dict(cds_targetinfo=cds_targetinfo, vi_class_input=vi_class_input,
                vi_class_labels=vi_class_labels, ifiberslider = ifiberslider,
                title=title, vi_file_fields = vi_file_fields, targetinfo=cds_targetinfo,
                model = cds_model, vi_z_input=vi_z_input, vi_category_select=vi_category_select),
        code=vi_class_code )
    vi_class_input.js_on_click(vi_class_callback)

    #- VI scanner name
    vi_name_input = TextInput(value=(cds_targetinfo.data['VI_scanner'][0]).strip(), title="Your name (3-letter acronym):")
    vi_name_code = js_files["CSVtoArray.js"] + js_files["save_vi.js"]
    vi_name_code += """
        for (var i=0; i<nspec; i++) {
            cds_targetinfo.data['VI_scanner'][i]=vi_name_input.value
        }
        var newname = vi_filename_input.value
        var pepe = newname.split("_")
        newname = ( pepe.slice(0,pepe.length-1).join("_") ) + ("_"+vi_name_input.value+".csv")
        vi_filename_input.value = newname
        autosave_vi_localStorage(vi_file_fields, cds_targetinfo.data, title)
        """
    vi_name_callback = CustomJS(
        args=dict(cds_targetinfo=cds_targetinfo, nspec = nspec, vi_name_input=vi_name_input,
                 vi_filename_input=vi_filename_input, title=title, vi_file_fields=vi_file_fields),
        code=vi_name_code )
    vi_name_input.js_on_change('value',vi_name_callback)

    #- Guidelines for VI flags
    vi_guideline_txt = "<B> VI guidelines </B>"
    vi_guideline_txt += "<BR /> <B> Classification flags: </B>"
    for flag in vi_flags :
        if flag['type'] == 'class' : vi_guideline_txt += ("<BR />&emsp;&emsp;[&emsp;"+flag['label']+"&emsp;] "+flag['description'])
    vi_guideline_txt += "<BR /> <B> Optional indications: </B>"
    for flag in vi_flags :
        if flag['type'] == 'issue' :
            vi_guideline_txt += ( "<BR />&emsp;&emsp;[&emsp;" + flag['label'] +
                                 "&emsp;(" + flag['shortlabel'] + ")&emsp;] " + flag['description'] )
    vi_guideline_txt += "<BR /> <B> Comments: </B> <BR /> 100 characters max, avoid commas (automatically replaced by semi-columns), ASCII only."
    vi_guideline_div = Div(text=vi_guideline_txt)

    #- Save VI info to CSV file
    save_vi_button = Button(label="Download VI", button_type="success")
    save_vi_code = js_files["FileSaver.js"] + js_files["CSVtoArray.js"] + js_files["save_vi.js"]
    save_vi_code += """
        download_vi_file(vi_file_fields, cds_targetinfo.data, vi_filename_input.value)
        """
    save_vi_callback = CustomJS(
        args=dict(cds_targetinfo=cds_targetinfo,
            vi_file_fields=vi_file_fields, vi_filename_input=vi_filename_input),
        code=save_vi_code )
    save_vi_button.js_on_event('button_click', save_vi_callback)

    #- Recover auto-saved VI data in browser
    recover_vi_button = Button(label="Recover auto-saved VI", button_type="default")
    recover_vi_code = js_files["CSVtoArray.js"] + js_files["recover_autosave_vi.js"]
    recover_vi_callback = CustomJS(
        args = dict(title=title, vi_file_fields=vi_file_fields, cds_targetinfo=cds_targetinfo,
                   ifiber=ifiberslider.value, vi_comment_input=vi_comment_input,
                   vi_name_input=vi_name_input, vi_class_input=vi_class_input, vi_issue_input=vi_issue_input,
                   vi_issue_slabels=vi_issue_slabels, vi_class_labels=vi_class_labels),
        code = recover_vi_code )
    recover_vi_button.js_on_event('button_click', recover_vi_callback)

    #- Clear all auto-saved VI
    clear_vi_button = Button(label="Clear all auto-saved VI", button_type="default")
    clear_vi_callback = CustomJS( args = dict(), code = """
        localStorage.clear()
        """ )
    clear_vi_button.js_on_event('button_click', clear_vi_callback)

    #- Show VI in a table
    vi_table_columns = [
        TableColumn(field="VI_class_flag", title="Flag", width=40),
        TableColumn(field="VI_issue_flag", title="Opt.", width=50),
        TableColumn(field="VI_z", title="VI z", width=50),
        TableColumn(field="VI_spectype", title="VI spectype", width=150),
        TableColumn(field="VI_comment", title="VI comment", width=200)
    ]
    vi_table = DataTable(source=cds_targetinfo, columns=vi_table_columns, width=500)
    vi_table.height = 10 * vi_table.row_height


    #-----
    #- Main js code to update plot
    #
    update_plot_code = (js_files["adapt_plotrange.js"] + js_files["interp_grid.js"] +
                        js_files["smooth_data.js"] + js_files["coadd_brz_cameras.js"] +
                        js_files["update_plot.js"])
    # ONGOING
    the_fit_results = None if template_dicts is None else template_dicts[1] # dirty
    update_plot = CustomJS(
        args = dict(
            spectra = cds_spectra,
            coaddcam_spec = cds_coaddcam_spec,
            model = cds_model,
            othermodel = cds_othermodel,
            model_2ndfit = cds_model_2ndfit,
            targetinfo = cds_targetinfo,
#            target_info_div = target_info_div,
## BYPASS DIV
            fit_results = the_fit_results,
            zcat_disp_cds = zcat_disp_cds,
            targ_disp_cds = targ_disp_cds,
            ifiberslider = ifiberslider,
            smootherslider = smootherslider,
            z_input = z_input,
            fig = fig,
            imfig_source=imfig_source,
            imfig_urls=imfig_urls,
            model_select = model_select,
            vi_comment_input = vi_comment_input,
            vi_std_comment_select = vi_std_comment_select,
            vi_name_input = vi_name_input,
            vi_class_input = vi_class_input,
            vi_class_labels = vi_class_labels,
            vi_issue_input = vi_issue_input,
            vi_z_input = vi_z_input,
            vi_category_select = vi_category_select,
            vi_issue_slabels = vi_issue_slabels
            ),
        code = update_plot_code
    )
    smootherslider.js_on_change('value', update_plot)
    ifiberslider.js_on_change('value', update_plot)


    #-----
    #- Bokeh setup
    # NB widget height / width are still partly hardcoded, but not arbitrary except for Spacers

    slider_width = plot_width - 2*navigation_button_width
    navigator = bl.row(
        bl.column(prev_button, width=navigation_button_width+15),
        bl.column(next_button, width=navigation_button_width+20),
        bl.column(ifiberslider, width=plot_width+(plot_height//2)-(60*len(vi_class_labels)+2*navigation_button_width+35))
    )
    if with_vi_widgets :
        navigator.children.insert(1, bl.column(vi_class_input, width=60*len(vi_class_labels)) )
        vi_widget_set = bl.column(
            bl.column( Div(text="VI optional indications :"), width=300 ),
            bl.row(
                bl.column(vi_issue_input, width=150),
                bl.column(vi_z_input, width=150),
                bl.column(vi_category_select, width=150)
            ),
            bl.row(
                bl.column(vi_comment_input, width=300),
                bl.column(vi_std_comment_select, width=200),
            ),
            bl.row(
                bl.column(vi_name_input, width=200),
                bl.column(vi_filename_input, width=300)
            ),
            bl.column(save_vi_button, width=100),
            bl.column(vi_table),
            bl.row(
                bl.column(recover_vi_button, width=150),
                bl.column(clear_vi_button, width=150)
            ),
            background='#f5f5f0'
        )
    plot_widget_set = bl.column(
        bl.column( Div(text="Pipeline fit: ") ),
        bl.column(zcat_display, width=plot_widget_width),
        bl.row(
            bl.column(
                bl.row(
                    bl.column(z_minus_button, width=z_button_width+15),
                    bl.column(zslider, width=plot_widget_width-2*z_button_width-135),
                    bl.column(z_plus_button, width=z_button_width)
                ),
                bl.row(
                    bl.column(dzslider, width=plot_widget_width-235),
                    bl.column(bl.Spacer(width=20)),
                    bl.column(zreset_button, width=100)
                )
            ),
            bl.column(bl.Spacer(width=15)),
            bl.column(
                bl.column(z_input, width=100),
                bl.column(z_tovi_button, width=100)
            ),
            background='#fff7e6'
        ),
        bl.column(smootherslider, width=plot_widget_width),
#        bl.column(display_options_group,width=120),
        bl.row(
            bl.column(coaddcam_buttons, width=200),
            bl.column(bl.Spacer(width=30)),
            bl.column(waveframe_buttons, width=120)
        ),
        bl.row(
            bl.column(lines_button_group, width=200),
            bl.column(bl.Spacer(width=30)),
            bl.column(majorline_checkbox, width=120)
        )
    )
    if model_select is not None :
        plot_widget_set.children.insert(3, bl.column(model_select, width=200))
    if with_vi_widgets :
        plot_widget_set.children.append( bl.column(bl.Spacer(height=30)) )
        plot_widget_set.children.append( bl.column(vi_guideline_div, width=plot_widget_width) )
        full_widget_set = bl.row(
            vi_widget_set,
            bl.column(bl.Spacer(width=40)),
            plot_widget_set
        )
    else : full_widget_set = plot_widget_set

    main_bokehsetup = bl.column(
        bl.row(fig, bl.column(imfig, zoomfig), bl.Spacer(width=20)),
        bl.row(
            bl.column(targ_display, width=600), # plot_width - 200
            bl.column(bl.Spacer(width=20)),
            bl.column(reset_plotrange_button, width = 120),
            bl.column(bl.Spacer(width=80)),
            bl.column(oii_zoom_button, width=80),
            bl.column(oii_undo_button, width=50),
        ),
        navigator,
        full_widget_set
    )

    if with_thumb_tab is False :
        full_viewer = main_bokehsetup
    else :
        full_viewer = Tabs()
        ncols_grid = 5 # TODO un-hardcode
        titles = None # TODO define
        miniplot_width = ( plot_width + (plot_height//2) ) // ncols_grid
        thumb_grid = grid_thumbs(spectra, miniplot_width, x_range=(xmin,xmax), ncols_grid=ncols_grid, titles=titles)
        tab1 = Panel(child = main_bokehsetup, title='Main viewer')
        tab2 = Panel(child = thumb_grid, title='Gallery')
        full_viewer.tabs=[ tab1, tab2 ]

        # Dirty trick : callback functions on thumbs need to be defined AFTER the full_viewer is implemented
        # Otherwise, at least one issue = no toolbar anymore for main fig. (apparently due to ifiberslider in callback args)
        for i_spec in range(nspec) :
            thumb_callback = CustomJS(args=dict(full_viewer=full_viewer, i_spec=i_spec, ifiberslider=ifiberslider), code="""
            full_viewer.active = 0
             ifiberslider.value = i_spec
            """)
            (thumb_grid.children[i_spec][0]).js_on_event(bokeh.events.DoubleTap, thumb_callback)

    if notebook:
        bk.show(full_viewer)
    else:
        bk.save(full_viewer)

    #-----
    #- "Light" Bokeh setup including only the thumbnail gallery
    if with_thumb_only_page :
        thumb_page = html_page.replace("specviewer_"+title, "thumbs_specviewer_"+title)
        bk.output_file(thumb_page, title='DESI spectral viewer - thumbnail gallery')
        ncols_grid = 5 # TODO un-hardcode
        titles = None # TODO define
        miniplot_width = ( plot_width + (plot_height//2) ) // ncols_grid
        thumb_grid = grid_thumbs(spectra, miniplot_width, x_range=(xmin,xmax), ncols_grid=ncols_grid, titles=titles)
        thumb_viewer = bl.column(
            bl.column( Div(text=
                           " <h3> Thumbnail gallery for DESI spectra in "+title+" </h3>" +
                           " <p> Click <a href='specviewer_"+title+".html'>here</a> to access the spectral viewer corresponding to these spectra. </p>"
                          ), width=plot_width ),
            bl.column( thumb_grid )
        )
        bk.save(thumb_viewer)

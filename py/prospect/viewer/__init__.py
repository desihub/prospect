# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
===============
prospect.viewer
===============

Run a spectral viewer (plot spectra and show widgets).

Spectra can be:

- DESI Spectra or Frames,
- `specutils`_-compatible objects (see :mod:`prospect.specutils` for objects and IO routines).

.. _`specutils`: https://specutils.readthedocs.io

"""

import os, sys
from pkg_resources import resource_filename

import numpy as np
import scipy.ndimage.filters

from astropy.table import Table
import astropy.io.fits

import bokeh.plotting as bk
from bokeh.models import ColumnDataSource, CDSView, IndexFilter
from bokeh.models import CustomJS, LabelSet, Label, Span, Legend, Panel, Tabs, BoxAnnotation
from bokeh.models.widgets import (
    Slider, Button, Div, CheckboxGroup, CheckboxButtonGroup, RadioButtonGroup,
    TextInput, Select, DataTable, TableColumn, Toggle)
import bokeh.layouts as bl

_specutils_imported = True
try:
    from specutils import Spectrum1D, SpectrumList
except ImportError:
    _specutils_imported = False

_desispec_imported = True
try:
    import desispec.io
    import desispec.spectra
    import desispec.frame
    from desispec.interpolation import resample_flux
except ImportError:
    _desispec_imported = False

_redrock_imported = True
try:
    import redrock.templates
    from redrock.archetypes import All_archetypes
except ImportError:
    _redrock_imported = False

from ..utilities import frames2spectra, create_zcat_from_redrock_cat
from .cds import ViewerCDS
from .plots import ViewerPlots
from .widgets import ViewerWidgets
from .vi_widgets import ViewerVIWidgets
from .layouts import ViewerLayout, StandaloneThumbLayout

def load_redrock_templates(template_dir=None) :
    '''
    Load redrock templates; redirect stdout because redrock is chatty
    '''
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
    return templates


def create_model(spectra, zbest, archetype_fit=False, archetypes_dir=None, template_dir=None):
    '''
    Returns model_wave[nwave], model_flux[nspec, nwave], row matched to zbest,
    which can be in a different order than spectra.
    - zbest must be entry-matched to spectra.
    '''

    assert _redrock_imported
    assert _desispec_imported  # for resample_flux

    if np.any(zbest['TARGETID'] != spectra.fibermap['TARGETID']) :
        raise ValueError('zcatalog and spectra do not match (different targetids)')

    if archetype_fit:
        archetypes = All_archetypes(archetypes_dir=archetypes_dir).archetypes
    else:
        templates = load_redrock_templates(template_dir=template_dir)

    #- Empty model flux arrays per band to fill
    model_flux = dict()
    for band in spectra.bands:
        model_flux[band] = np.zeros(spectra.flux[band].shape)

    for i in range(len(zbest)):
        zb = zbest[i]

        if archetype_fit:
            archetype  = archetypes[zb['SPECTYPE']]
            coeff      = zb['COEFF']
            
            for band in spectra.bands:
                wave                = spectra.wave[band]
                wavehash            = hash((len(wave), wave[0], wave[1], wave[-2], wave[-1], spectra.R[band].data.shape[0]))
                dwave               = {wavehash: wave}
                mx                  = archetype.eval(zb['SUBTYPE'], dwave, coeff, wave, zb['Z']) * (1+zb['Z'])
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
        raise RuntimeError("create_model: Set of bands for spectra not supported")

    return model_wave, mflux



def make_template_dicts(redrock_cat, delta_lambd_templates=3, with_fit_templates=True, template_dir=None) :
    """
    Input : TODO document
    - redrock_cat : Table produced by match_redrock_zfit_to_spectra (matches spectra).
    Create list of CDS including all data needed to plot other models (Nth best fit, std templates) :
    - list of templates used in fits
    - RR output for Nth best fits
    - list of std templates
    """

    assert _redrock_imported
    assert _desispec_imported # for resample_flux
    rr_templts = load_redrock_templates(template_dir=template_dir)

    if with_fit_templates :
        dict_fit_templates = dict()
        for key,val in rr_templts.items() :
            fulltype_key = "_".join(key)
            wave_array = np.arange(val.wave[0],val.wave[-1],delta_lambd_templates)
            flux_array = np.zeros(( val.flux.shape[0],len(wave_array) ))
            for i in range(val.flux.shape[0]) :
                flux_array[i,:] = resample_flux(wave_array, val.wave, val.flux[i,:])
            dict_fit_templates["wave_"+fulltype_key] = wave_array
            dict_fit_templates["flux_"+fulltype_key] = flux_array
    else : dict_fit_templates = None

    dict_fit_results = dict()
    for key in redrock_cat.keys() :
        dict_fit_results[key] = np.asarray(redrock_cat[key])
    dict_fit_results['Nfit'] = redrock_cat['Z'].shape[1]

    # TODO fix the list of std templates
    # We take flux[0,:] : ie use first entry in RR template basis
    # We choose here not to convolve with a "typical" resolution (could easily be done)
    # Std template : corresponding RR template . TODO put this list somewhere else
    std_templates = {'QSO': ('QSO',''), 'GALAXY': ('GALAXY',''), 'STAR': ('STAR','F') }
    dict_std_templates = dict()
    for key,rr_key in std_templates.items() :
        wave_array = np.arange(rr_templts[rr_key].wave[0],rr_templts[rr_key].wave[-1],delta_lambd_templates)
        flux_array = resample_flux(wave_array, rr_templts[rr_key].wave, rr_templts[rr_key].flux[0,:])
        dict_std_templates["wave_"+key] = wave_array
        dict_std_templates["flux_"+key] = flux_array

    return [dict_fit_templates, dict_fit_results, dict_std_templates]


def plotspectra(spectra, zcatalog=None, redrock_cat=None, notebook=False, html_dir=None, title=None,
                with_imaging=True, with_noise=True, with_thumb_tab=True, with_vi_widgets=True,
                top_metadata=None, vi_countdown=-1, with_thumb_only_page=False,
                with_coaddcam=True, mask_type='DESI_TARGET',
                model_from_zcat=True, model=None, num_approx_fits=None, with_full_2ndfit=True,
                template_dir=None, archetype_fit=False, archetypes_dir=None):
    '''Main prospect routine. From a set of spectra, creates a bokeh document
    used for VI, to be displayed as an HTML page or within a Jupyter notebook.

    Parameters
    ----------
    spectra : :class:`~desispec.spectra.Spectra` or :class:`~specutils.Spectrum1D` or :class:`~specutils.SpectrumList` or list of :class:`~desispec.frame.Frame`
        Input spectra. :class:`~specutils.Spectrum1D` are assumed to be SDSS/BOSS/eBOSS.
        Otherwise DESI spectra or frames is assumed.
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
        If ``False``, don't include thumb image from https://www.legacysurvey.org/viewer.
    with_noise : :class:`bool`, optional
        If ``False``, don't include uncertainty for each spectrum.
    with_thumb_tab : :class:`bool`, optional
        If ``False``, don't include a tab with spectra thumbnails.
    with_vi_widgets : :class:`bool`, optional
        Include widgets used to enter VI information. Set it to ``False`` if
        you do not intend to record VI files.
    top_metadata : :class:`list`, optional
        List of metadata to be highlighted in the top (most visible) table.
        Default values ['TARGETID', 'EXPID']
    vi_countdown : :class:`int`, optional
        If ``>0``, add a countdown widget in the VI panel, with a value in minutes given
        by `vi_countdown``.
    with_thumb_only_page : :class:`bool`, optional
        When creating a static HTML (`notebook` is ``False``), a light HTML
        page including only the thumb gallery will also be produced.
    with_coaddcam : :class:`bool`, optional
        Include camera-coaddition, only relevant for DESI.
    mask_type : :class:`str`, optional (default: DESI_TARGET)
        Bitmask type to identify target categories in the spectra. For DESI
        these could be: DESI_TARGET, SV1_DESI_TARGET, SV1_BGS_TARGET, CMX_TARGET.
    model_from_zcat : :class:`bool`, optional
        If ``True``, model spectra will be computed from the input `zcatalog`.
    model : :func:`tuple`, optional
        If set, use this input set of model spectra instead of computing it from `zcatalog`.
        model consists of (mwave, mflux); model must be entry-matched to `zcatalog`.
    num_approx_fits : :class:`int`, optional
        Number of best-fit models to display, if `redrock_cat` is provided.
        By default, all best-fit models available in `redrock_cat` are diplayed.
    with_full_2ndfit : :class:`bool`, optional
        If ``True``, the second best-fit model from `redrock_cat` will be displayed
        without approximation (no undersampling, full resolution).
    template_dir : :class:`str`, optional
        Redrock template directory.
    archetype_fit : :class:`bool`, optional
        If ``True``, assume `zcatalog` derived from :command:`redrock --archetypes`
        and plot model accordingly.
    archetypes_dir : :class:`str`, optional
        Directory path for archetypes if not :envvar:`RR_ARCHETYPE_DIR`.
    '''

    #- Check input spectra.
    #- Set masked bins to NaN for compatibility with bokeh.
    if _specutils_imported and isinstance(spectra, Spectrum1D):
        # We will assume this is from an SDSS/BOSS/eBOSS spPlate file.
        survey = 'SDSS'
        nspec = spectra.flux.shape[0]
        bad = (spectra.uncertainty.array == 0.0) | spectra.mask
        spectra.flux[bad] = np.nan
    elif _specutils_imported and isinstance(spectra, SpectrumList):
        # We will assume this is from a DESI spectra-64 file.
        survey = 'DESI'
        nspec = spectra[0].flux.shape[0]
        for s in spectra:
            bad = (s.uncertainty.array == 0.0) | s.mask
            s.flux[bad] = np.nan
    else:
        # DESI object (Spectra or list of Frame)
        survey = 'DESI'
        if _desispec_imported and isinstance(spectra, desispec.spectra.Spectra):
            nspec = spectra.num_spectra()
        elif _desispec_imported and isinstance(spectra, list) and isinstance(spectra[0], desispec.frame.Frame):
        # If inputs are frames, convert to a spectra object
            spectra = frames2spectra(spectra)
            nspec = spectra.num_spectra()
            if title is None:
                title = 'Night {} ExpID {} Spectrograph {}'.format(
                    spectra.meta['NIGHT'], spectra.meta['EXPID'], spectra.meta['CAMERA'][1],
                )
        else:
            raise ValueError("Unsupported type for input spectra. \n"+
                    "    _specutils_imported = "+str(_specutils_imported)+"\n"+ 
                    "    _desispec_imported = "+str(_desispec_imported))
        for band in spectra.bands:
            bad = (spectra.ivar[band] == 0.0) | (spectra.mask[band] != 0)
            spectra.flux[band][bad] = np.nan
        #- No coaddition if spectra is already single-band
        if len(spectra.bands)==1 : with_coaddcam = False
    
    if title is None:
        title = "specviewer"

    #- Input zcatalog / model
    if zcatalog is not None:
        if survey == 'SDSS':
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

        if model_from_zcat :
            # DESI spectra will obtain the model from templates.
            model = create_model(spectra, zcatalog,
                                 archetype_fit=archetype_fit,
                                 archetypes_dir=archetypes_dir,
                                 template_dir=template_dir)

    #-----
    #- Gather information into ColumnDataSource objects for Bokeh
    viewer_cds = ViewerCDS()
    viewer_cds.load_spectra(spectra, with_noise)
    if with_coaddcam :
        viewer_cds.init_coaddcam_spec(spectra, with_noise)
    if model is not None:
        viewer_cds.init_model(model)

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
            viewer_cds.init_model(model_2ndfit, second_fit=True)
        viewer_cds.init_othermodel(zcatalog)
    else :
        template_dicts = None

    viewer_cds.load_metadata(spectra, mask_type=mask_type, zcatalog=zcatalog, survey=survey)
    
    #-------------------------
    #-- Graphical objects --
    #-------------------------

    viewer_plots = ViewerPlots()
    viewer_plots.create_mainfig(spectra, title, viewer_cds, survey,
                                with_noise=with_noise, with_coaddcam=with_coaddcam)
    viewer_plots.create_zoomfig(viewer_cds, 
                                with_noise=with_noise, with_coaddcam=with_coaddcam)
    if with_imaging :
        viewer_plots.create_imfig(spectra)

    #-----
    #- Emission and absorption lines
    z = zcatalog['Z'][0] if (zcatalog is not None) else 0.0
    viewer_cds.load_spectral_lines(z)
    viewer_plots.add_spectral_lines(viewer_cds, figure='main')
    viewer_plots.add_spectral_lines(viewer_cds, figure='zoom', label_offsets=[50, 5])
    

    #-------------------------
    #-- Widgets and callbacks --
    #-------------------------

    viewer_widgets = ViewerWidgets(viewer_plots, nspec)
    viewer_widgets.add_navigation(nspec)
    viewer_widgets.add_resetrange(viewer_cds, viewer_plots)

    viewer_widgets.add_redshift_widgets(z, viewer_cds, viewer_plots)
    viewer_widgets.add_oii_widgets(viewer_plots)

    viewer_plots.add_imfig_callback(viewer_widgets)

    if viewer_cds.cds_coaddcam_spec is not None :
        viewer_widgets.add_coaddcam(viewer_plots)

    if zcatalog is not None :
        show_zcat = True
    else : show_zcat = False
    if top_metadata is None: top_metadata = ['TARGETID', 'EXPID']
    viewer_widgets.add_metadata_tables(viewer_cds, top_metadata=top_metadata,
                                       show_zcat=show_zcat, template_dicts=template_dicts)
    viewer_widgets.add_specline_toggles(viewer_cds, viewer_plots)

    if template_dicts is not None :
        viewer_cds.compute_median_spectra(spectra)
        viewer_widgets.add_model_select(viewer_cds, template_dicts, 
                            num_approx_fits, with_full_2ndfit=with_full_2ndfit)
    
    #-----
    #- VI-related widgets
    ## TODO if with_vi_widgets (need to adapt update_plot.js..)

    viewer_vi_widgets = ViewerVIWidgets(title, viewer_cds)

    viewer_vi_widgets.add_filename()
    viewer_vi_widgets.add_vi_issues(viewer_cds, viewer_widgets)
    viewer_vi_widgets.add_vi_z(viewer_cds, viewer_widgets)
    viewer_vi_widgets.add_vi_spectype(viewer_cds, viewer_widgets)
    viewer_vi_widgets.add_vi_comment(viewer_cds, viewer_widgets)
    viewer_vi_widgets.add_vi_classification(viewer_cds, viewer_widgets)
    viewer_vi_widgets.add_vi_scanner(viewer_cds, nspec)
    viewer_vi_widgets.add_guidelines()
    viewer_vi_widgets.add_vi_storage(viewer_cds, viewer_widgets)
    viewer_vi_widgets.add_vi_table(viewer_cds)

    if (vi_countdown > 0) :
        viewer_vi_widgets.add_countdown(vi_countdown)

    viewer_widgets.add_update_plot_callback(viewer_cds, viewer_plots, 
                viewer_vi_widgets, template_dicts)

    #-----
    #- Bokeh layout and output
    
    bokeh_layout = ViewerLayout(viewer_plots, viewer_widgets, viewer_vi_widgets,
                              with_vi_widgets=with_vi_widgets)
    if with_thumb_tab:
        bokeh_layout.add_thumb_tab(spectra, viewer_plots, viewer_widgets, nspec)

    if notebook:
        bk.output_notebook()
        bk.show(bokeh_layout.full_viewer)
    else:
        if html_dir is None : raise RuntimeError("Need html_dir")
        html_page = os.path.join(html_dir, "specviewer_"+title+".html")
        bk.output_file(html_page, title='DESI spectral viewer')
        bk.save(bokeh_layout.full_viewer)

    #-----
    #- "Light" Bokeh layout including only the thumbnail gallery
    if with_thumb_only_page :
        assert not notebook
        thumb_page = os.path.join(html_dir, "thumbs_specviewer_"+title+".html")
        bk.output_file(thumb_page, title='DESI spectral viewer - thumbnail gallery')
        thumb_grid = StandaloneThumbLayout(spectra, viewer_plots, title)
        bk.save(thumb_grid.thumb_viewer)

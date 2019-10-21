# -*- coding: utf-8 -*-


"""
TODO
* add target details tab
* add code details tab (version, SPECPROD)
* redshift model fit
* better smoothing kernel, e.g. gaussian
* sky? ivar? (could be useful, but makes data payload larger)
"""

import os, sys
import argparse

import numpy as np
from astropy.table import Table
import astropy.io.fits

import bokeh.plotting as bk
from bokeh.models import ColumnDataSource, CDSView, IndexFilter
from bokeh.models import CustomJS, LabelSet, Label, Span, Legend
from bokeh.models.widgets import (
    Slider, Button, Div, CheckboxGroup, CheckboxButtonGroup, RadioButtonGroup, TextInput, Select)
from bokeh.layouts import widgetbox
import bokeh.events
# from bokeh.layouts import row, column

import desispec.io
from desitarget.targetmask import desi_mask
from desitarget.sv1.sv1_targetmask import desi_mask as desi_mask_sv1
import desispec.spectra
import desispec.frame

#from . import utils_specviewer
from prospect import utils_specviewer
from prospect import mycoaddcam
from astropy.table import Table

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


def frames2spectra(frames):
    '''Convert input list of Frames into Spectra object

    Do no propagate resolution, scores
    '''
    bands = list()
    wave = dict()
    flux = dict()
    ivar = dict()
    mask = dict()
    for fr in frames:
        band = fr.meta['CAMERA'][0]
        bands.append(band)
        wave[band] = fr.wave
        flux[band] = fr.flux
        ivar[band] = fr.ivar
        mask[band] = fr.mask

    spectra = desispec.spectra.Spectra(
        bands, wave, flux, ivar, mask, fibermap=fr.fibermap, meta=fr.meta
    )
    return spectra

def create_model(spectra, zbest):
    '''
    Returns model_wave[nwave], model_flux[nspec, nwave], row matched to zbest,
    which can be in a different order than spectra.
    NB currently, zbest must have the same size as spectra.
    '''
    import redrock.templates
    from desispec.interpolation import resample_flux

    nspec = spectra.num_spectra()
    assert len(zbest) == nspec

    #- Load redrock templates; redirect stdout because redrock is chatty
    saved_stdout = sys.stdout
    sys.stdout = open('/dev/null', 'w')
    try:
        templates = dict()
        for filename in redrock.templates.find_templates():
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

    targetids = spectra.target_ids()
    for i in range(len(zbest)):
        zb = zbest[i]
        j = np.where(targetids == zb['TARGETID'])[0][0]

        tx = templates[(zb['SPECTYPE'], zb['SUBTYPE'])]
        coeff = zb['COEFF'][0:tx.nbasis]
        model = tx.flux.T.dot(coeff).T
        for band in spectra.bands:
            mx = resample_flux(spectra.wave[band], tx.wave*(1+zb['Z']), model)
            model_flux[band][i] = spectra.R[band][j].dot(mx)

    #- Now combine to a single wavelength grid across all cameras
    #- TODO: assumes b,r,z all exist
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

    return model_wave, mflux


def _viewer_urls(spectra, zoom=13, layer='dr8'):
    """Return legacysurvey.org viewer URLs for all spectra.
    """
    u = "http://legacysurvey.org/viewer/jpeg-cutout?ra={0:f}&dec={1:f}&zoom={2:d}&layer={3}"
    v = "http://legacysurvey.org/viewer/?ra={0:f}&dec={1:f}&zoom={2:d}&layer={3}"
    try:
        ra = spectra.fibermap['RA_TARGET']
        dec = spectra.fibermap['DEC_TARGET']
    except KeyError:
        ra = spectra.fibermap['TARGET_RA']
        dec = spectra.fibermap['TARGET_DEC']

    return [(u.format(ra[i], dec[i], zoom, layer),
             v.format(ra[i], dec[i], zoom, layer),
             'RA, Dec = {0:.4f}, {1:+.4f}'.format(ra[i], dec[i]))
            for i in range(len(ra))]


def plotspectra(spectra, zcatalog=None, model_from_zcat=True, model=None, notebook=False, vidata=None, savedir='.', is_coadded=True, title=None, html_dir=None, with_noise=True, sv=False):
    '''
    Main prospect routine, creates a bokeh document from a set of spectra and fits

    Parameter
    ---------
    spectra : desi spectra object, or a list of frames
    zcatalog : FITS file of pipeline redshifts for the spectra. Currently supports only redrock-PCA files.
    model_from_zcat : if True, model spectra will be computed from the input zcatalog
    model : if set, use this input set of model spectra (instead of computing it from zcat)
        model format (mwave, mflux); model must be entry-matched to zcatalog.
    notebook : if True, bokeh outputs the viewer to notebook, else to a (static) html page
    vidata : VI information to be preloaded and displayed. Currently disabled.
    is_coadded : set to True if spectra are coadds
    title : title of produced html page and bokeh figure
    html_dir : directory to store html page
    with_noise : include noise for each spectrum
    sv : if True, will use SV1_DESI_TARGET instead of DESI_TARGET
    '''

    #- If inputs are frames, convert to a spectra object
    if isinstance(spectra, list) and isinstance(spectra[0], desispec.frame.Frame):
        spectra = frames2spectra(spectra)
        frame_input = True
    else:
        frame_input = False

    if frame_input and title is None:
        meta = spectra.meta
        title = 'Night {} ExpID {} Spectrograph {}'.format(
            meta['NIGHT'], meta['EXPID'], meta['CAMERA'][1],
        )

    if notebook:
        bk.output_notebook()
    else :
        if title is None : title="specviewer"
        if html_dir is None : raise RuntimeError("Need html_dir")
        bk.output_file(html_dir+"/"+title+".html", title='DESI spectral viewer')

    #- Gather spectra into ColumnDataSource objects for Bokeh
    nspec = spectra.num_spectra()
    cds_spectra = list()

    for band in spectra.bands:
        #- Set masked bins to NaN so that Bokeh won't plot them
        bad = (spectra.ivar[band] == 0.0) | (spectra.mask[band] != 0)
        spectra.flux[band][bad] = np.nan

        cdsdata=dict(
            origwave=spectra.wave[band].copy(),
            plotwave=spectra.wave[band].copy(),
            )

        for i in range(nspec):
            key = 'origflux'+str(i)
            cdsdata[key] = spectra.flux[band][i]
            if with_noise :
                key = 'orignoise'+str(i)
                noise = np.zeros(len(spectra.ivar[band][i]))
                w, = np.where( (spectra.ivar[band][i] > 0))
                noise[w] = 1/np.sqrt(spectra.ivar[band][i][w])
                cdsdata[key] = noise

        cdsdata['plotflux'] = cdsdata['origflux0']
        if with_noise : cdsdata['plotnoise'] = cdsdata['orignoise0'] 

        cds_spectra.append(
            bk.ColumnDataSource(cdsdata, name=band)
            )
    
    # CDS object for camera-merged spectrum
    # Do NOT store all coadded spectra in CDS obj, to reduce size of html files
    # Except for the first spectrum, coaddition is done later in javascript
    coadd_wave, coadd_flux, coadd_ivar = mycoaddcam.mycoaddcam(spectra)
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

    #- Reorder zcatalog to match input targets
    #- TODO: allow more than one zcatalog entry with different ZNUM per targetid
    if zcatalog is not None:
        ## (at the moment, keep argsorts-based code for record)
        #targetids = spectra.target_ids()
        #ii = np.argsort(np.argsort(targetids))
        #jj = np.argsort(zcatalog['TARGETID'])
        #kk = jj[ii]
        #zcatalog = zcatalog[kk]
        ##- That sequence of argsorts may feel like magic,
        ##- so make sure we got it right
        #assert np.all(zcatalog['TARGETID'] == targetids)
        #assert np.all(zcatalog['TARGETID'] == spectra.fibermap['TARGETID'])
        zcatalog, kk = utils_specviewer.match_zcat_to_spectra(zcatalog, spectra)
        
        #- Also need to re-order input model fluxes
        if model is not None :
            assert(model_from_zcat == False)
            mwave, mflux = model
            model = mwave, mflux[kk]

        if model_from_zcat == True :
            model = create_model(spectra, zcatalog)
            
    #- Gather models into ColumnDataSource objects, row matched to spectra
    if model is not None:
        mwave, mflux = model
        model_obswave = mwave.copy()
        model_restwave = mwave.copy()
        cds_model_data = dict(
            origwave = mwave.copy(),
            plotwave = mwave.copy(),
            plotflux = np.zeros(len(mwave)),
        )

        for i in range(nspec):
            key = 'origflux'+str(i)
            cds_model_data[key] = mflux[i]

        cds_model_data['plotflux'] = cds_model_data['origflux0']
        cds_model = bk.ColumnDataSource(cds_model_data)
    else:
        cds_model = None

    #- Subset of zcatalog and fibermap columns into ColumnDataSource
    target_info = list()
    vi_info = list()
    for i, row in enumerate(spectra.fibermap):
        if sv :
            target_bit_names = ' '.join(desi_mask_sv1.names(row['SV1_DESI_TARGET']))
        else :
            target_bit_names = ' '.join(desi_mask.names(row['DESI_TARGET']))
        txt = 'Target {}: {}'.format(row['TARGETID'], target_bit_names)
        if zcatalog is not None:
            txt += '<BR/>{} z={:.4f} ± {:.4f}  ZWARN={}'.format(
                zcatalog['SPECTYPE'][i],
                zcatalog['Z'][i],
                zcatalog['ZERR'][i],
                zcatalog['ZWARN'][i],
            )
        target_info.append(txt)
        if ( (vidata is not None) and (len(vidata[i])>0) ) :
            txt_viinfo = ('<BR/> VI info : SCANNER FLAG COMMENTS')
            for the_vi in vidata[i] :
                txt_viinfo += ('<BR/>&emsp;&emsp;&emsp;&emsp; {0} {1} {2}'.format(the_vi['scannername'], the_vi['scanflag'], the_vi['VIcomment']))
        else : txt_viinfo = ('<BR/> No VI previously recorded for this target')
        vi_info.append(txt_viinfo)

    cds_targetinfo = bk.ColumnDataSource(
        dict(target_info=target_info),
        name='targetinfo')
    cds_targetinfo.add(vi_info, name='vi_info')
    if zcatalog is not None:
        cds_targetinfo.add(zcatalog['Z'], name='z')
        cds_targetinfo.add(zcatalog['SPECTYPE'].astype('U{0:d}'.format(zcatalog['SPECTYPE'].dtype.itemsize)), name='spectype')
    username = '-'
    if notebook and ("USER" in os.environ) : username = os.environ['USER']
    cds_targetinfo.add([username for i in range(nspec)], name='VI_ongoing_scanner')
    cds_targetinfo.add(['-' for i in range(nspec)], name='VI_ongoing_flag')
    cds_targetinfo.add(['-' for i in range(nspec)], name='VI_ongoing_comment')
    if not is_coadded :
        cds_targetinfo.add(spectra.fibermap['EXPID'], name='expid')
        cds_targetinfo.add(spectra.fibermap['FIBER'], name='fiber')
    else : # If coadd, fill VI accordingly
        cds_targetinfo.add(['-1' for i in range(nspec)], name='expid')
        cds_targetinfo.add(['-1' for i in range(nspec)], name='fiber')
    cds_targetinfo.add([str(x) for x in spectra.fibermap['TARGETID']], name='targetid') # !! No int64 in js !!

    #- FIXME: should not hardcode which DEPVERnn has which versions
    ### cds_targetinfo.add([spectra.meta['DEPVER10'] for i in range(nspec)], name='spec_version')
    ### cds_targetinfo.add([spectra.meta['DEPVER13'] for i in range(nspec)], name='redrock_version')
    cds_targetinfo.add(np.zeros(nspec), name='spec_version')
    cds_targetinfo.add(np.zeros(nspec), name='redrock_version')

    #- Determine initial ymin, ymax, xmin, xmax
    ymin = ymax = xmax = 0.0
    xmin = 100000.
    xmargin = 300.
    for band in spectra.bands:
        ymin = min(ymin, np.min(spectra.flux[band][0]))
        ymax = max(ymax, np.max(spectra.flux[band][0]))
        xmin = min(xmin, np.min(spectra.wave[band]))
        xmax = max(xmax, np.max(spectra.wave[band]))        
    xmin -= xmargin
    xmax += xmargin
        
    plot_width=800
    plot_height=400
    # tools = 'pan,box_zoom,wheel_zoom,undo,redo,reset,save'
    tools = 'pan,box_zoom,wheel_zoom,save' # reset
    fig = bk.figure(height=plot_height, width=plot_width, title=title,
        tools=tools, toolbar_location='above', y_range=(ymin, ymax), x_range=(xmin, xmax))
    fig.toolbar.active_drag = fig.tools[0]    #- pan zoom (previously box)
    fig.toolbar.active_scroll = fig.tools[2]  #- wheel zoom
    fig.xaxis.axis_label = 'Wavelength [Å]'
    fig.yaxis.axis_label = 'Flux'
    fig.xaxis.axis_label_text_font_style = 'normal'
    fig.yaxis.axis_label_text_font_style = 'normal'
    colors = dict(b='#1f77b4', r='#d62728', z='maroon', coadd='#d62728')
    noise_colors = dict(b='greenyellow', r='green', z='forestgreen', coadd='green') # TODO test several and choose
    alpha_discrete = 0.2 # alpha for "almost-hidden" curves (single-arm spectra and noise by default)
    
    data_lines = list()
    for spec in cds_spectra:
        lx = fig.line('plotwave', 'plotflux', source=spec, line_color=colors[spec.name], line_alpha=alpha_discrete)
        data_lines.append(lx)
    lx = fig.line('plotwave', 'plotflux', source=cds_coaddcam_spec, line_color=colors['coadd'], line_alpha=1)
    data_lines.append(lx)
    
    noise_lines = list()
    if with_noise :
        for spec in cds_spectra :
            lx = fig.line('plotwave', 'plotnoise', source=spec, line_color=noise_colors[spec.name], line_alpha=alpha_discrete)
            noise_lines.append(lx)
        lx = fig.line('plotwave', 'plotnoise', source=cds_coaddcam_spec, line_color=noise_colors['coadd'], line_alpha=1)
        noise_lines.append(lx)

    model_lines = list()
    if cds_model is not None:
        lx = fig.line('plotwave', 'plotflux', source=cds_model, line_color='black')
        model_lines.append(lx)

    legend_items = [("data",  data_lines[-1::-1])] #- reversed to get blue as lengend entry
    if cds_model is not None : 
        legend_items.append(("model", model_lines))
    if with_noise : 
        legend_items.append(("noise", noise_lines[-1::-1])) # same as for data_lines
    legend = Legend(items=legend_items)

    fig.add_layout(legend, 'center')
    fig.legend.click_policy = 'hide'    #- or 'mute'

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
    zoom_data_lines.append(zoomfig.line('plotwave', 'plotflux', source=cds_coaddcam_spec, line_color=colors['coadd'], line_alpha=1))
    
    if with_noise :
        lx = zoomfig.line('plotwave', 'plotnoise', source=cds_coaddcam_spec, line_color=noise_colors['coadd'], line_alpha=1)
        zoom_noise_lines.append(lx)
            
    zoom_model_lines = list()
    if cds_model is not None:
        zoom_model_lines.append(zoomfig.line('plotwave', 'plotflux', source=cds_model, line_color='black'))

    #- Callback to update zoom window x-range
    zoom_callback = CustomJS(
        args=dict(zoomfig=zoomfig,fig=fig),
        code="""
            zoomfig.x_range.start = cb_obj.x - 100;
            zoomfig.x_range.end = cb_obj.x + 100;
        """)

    fig.js_on_event(bokeh.events.MouseMove, zoom_callback)

    #
    # Targeting image
    #
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

    #-----
    #- Emission and absorption lines
    z = zcatalog['Z'][0] if (zcatalog is not None) else 0.0
    line_data, lines, line_labels = add_lines(fig, z=z)
    zoom_line_data, zoom_lines, zoom_line_labels = add_lines(zoomfig, z=z, label_offsets=[50, 5])

    #-----
    #- Add widgets for controling plots
    z1 = np.floor(z*100)/100
    dz = z-z1
    zslider = Slider(start=0.0, end=4.0, value=z1, step=0.01, title='Redshift')
    dzslider = Slider(start=-0.01, end=0.01, value=dz, step=0.0001, title='+ Delta redshift')
    dzslider.format = "0[.]0000"
    z_display = Div(text="<br><b>z<sub>disp</sub> = "+("{:.4f}").format(z+dz)+"</b>")

    #- Observer vs. Rest frame wavelengths
    waveframe_buttons = RadioButtonGroup(
        labels=["Obs", "Rest"], active=0)

    ifiberslider = Slider(start=0, end=nspec-1, value=0, step=1)
    ifiberslider.title = 'Spectrum'
#     if frame_input:
#         ifiberslider.title = 'Fiber'
#     else:
#         ifiberslider.title = 'Target'

    zslider_callback  = CustomJS(
        args=dict(
            spectra = cds_spectra,
            coaddcam_spec = cds_coaddcam_spec,
            model = cds_model,
            targetinfo = cds_targetinfo,
            ifiberslider = ifiberslider,
            zslider=zslider,
            dzslider=dzslider,
            z_display = z_display,
            waveframe_buttons=waveframe_buttons,
            line_data=line_data, lines=lines, line_labels=line_labels,
            zlines=zoom_lines, zline_labels=zoom_line_labels,
            fig=fig,
            ),
        code="""
        var z = zslider.value + dzslider.value
        z_display.text = "<br><b>z<sub>disp</sub> = " + z.toFixed(4) + "</b>"
        var line_restwave = line_data.data['restwave']
        var ifiber = ifiberslider.value
        var zfit = 0.0
        if(targetinfo.data['z'] != undefined) {
            zfit = targetinfo.data['z'][ifiber]
        }
        for(var i=0; i<line_restwave.length; i++) {
            var waveshift = (waveframe_buttons.active == 0) ? 1+z : 1 ;
            lines[i].location = line_restwave[i] * waveshift ;
            line_labels[i].x = line_restwave[i] * waveshift ;
            zlines[i].location = line_restwave[i] * waveshift ;
            zline_labels[i].x = line_restwave[i] * waveshift ;
        }
        
        function shift_plotwave(cds_spec, waveshift) {
            var data = cds_spec.data
            var origwave = data['origwave']
            var plotwave = data['plotwave']
            for (var j=0; j<plotwave.length; j++) {
                plotwave[j] = origwave[j] * waveshift ;
            }
            cds_spec.change.emit()
        }
        
        var waveshift_spec = (waveframe_buttons.active == 0) ? 1 : 1/(1+z) ;
        for(var i=0; i<spectra.length; i++) {
            shift_plotwave(spectra[i], waveshift_spec)
        }
        shift_plotwave(coaddcam_spec, waveshift_spec)
        
        // Update model wavelength array
        if(model) {
            var waveshift_model = (waveframe_buttons.active == 0) ? (1+z)/(1+zfit) : 1/(1+zfit) ;
            shift_plotwave(model, waveshift_model)
        }
        """)

    zslider.js_on_change('value', zslider_callback)
    dzslider.js_on_change('value', zslider_callback)
    waveframe_buttons.js_on_click(zslider_callback)

    zreset_button = Button(label='Reset redshift')
    zreset_callback = CustomJS(
        args=dict(zslider=zslider, dzslider=dzslider, targetinfo=cds_targetinfo, ifiberslider=ifiberslider),
        code="""
            var ifiber = ifiberslider.value
            var z = targetinfo.data['z'][ifiber]
            var z1 = Math.floor(z*100) / 100
            zslider.value = z1
            dzslider.value = (z - z1)
        """)
    zreset_button.js_on_event('button_click', zreset_callback)

    plotrange_callback = CustomJS(
        args = dict(
            zslider=zslider,
            dzslider=dzslider,
            waveframe_buttons=waveframe_buttons,
            fig=fig,
        ),
        code="""
        var z = zslider.value + dzslider.value
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

    imfig_callback = CustomJS(args=dict(urls=imfig_urls,
                                        ifiberslider=ifiberslider),
                              code='''window.open(urls[ifiberslider.value][1], "_blank");''')
    imfig.js_on_event('tap', imfig_callback)

    smootherslider = Slider(start=0, end=31, value=0, step=1.0, title='Gaussian Sigma Smooth')
   
    #-----
    #- Checkboxes to display or not noise and model
    if with_noise : 
        display_options_group = CheckboxGroup(labels=['Show model', 'Show noise spectra'], active=[0,1])
    else :
        display_options_group = CheckboxGroup(labels=['Show model'], active=[0])
    disp_opt_callback = CustomJS(
        args = dict(noise_lines=noise_lines, model_lines=model_lines, zoom_noise_lines=zoom_noise_lines, zoom_model_lines=zoom_model_lines), code="""
        for (var i=0; i<noise_lines.length; i++) {
            if (cb_obj.active.indexOf(1) >= 0) {
                noise_lines[i].visible = true
                zoom_noise_lines[i].visible = true
            } else {
                noise_lines[i].visible = false
                zoom_noise_lines[i].visible = false
            }
        }
        for (var i=0; i<model_lines.length; i++) {
            if (cb_obj.active.indexOf(0) >= 0) {
                model_lines[i].visible = true
                zoom_model_lines[i].visible = true
            } else {
                model_lines[i].visible = false
                zoom_model_lines[i].visible = false
            }
        }
        """
    )
    display_options_group.js_on_click(disp_opt_callback)
    
    #-----
    #- Highlight individual-arm or camera-coadded spectra
    coaddcam_buttons = RadioButtonGroup( labels=["Camera-coadded spectrum", "Single-arm spectra"], active=0 )
    coaddcam_callback = CustomJS(
        args = dict(coaddcam_buttons=coaddcam_buttons, list_lines=[data_lines, noise_lines, zoom_data_lines, zoom_noise_lines], alpha_discrete=alpha_discrete), code="""
        var n_lines = list_lines[0].length
        for (var i=0; i<n_lines; i++) {
            var new_alpha = 1
            if (coaddcam_buttons.active == 0 && i<n_lines-1) new_alpha = alpha_discrete
            if (coaddcam_buttons.active == 1 && i==n_lines-1) new_alpha = alpha_discrete
            for (var j=0; j<list_lines.length; j++) {
                list_lines[j][i].glyph.line_alpha = new_alpha
            }
        }
        """
    )
    coaddcam_buttons.js_on_click(coaddcam_callback)

    target_info_div = Div(text=target_info[0])
    vi_info_div = Div(text=" ") # consistent with show_prev_vi="No" by default

    #-----
    #- Toggle lines
    lines_button_group = CheckboxButtonGroup(
            labels=["Emission", "Absorption"], active=[])

    lines_callback = CustomJS(
        args = dict(line_data=line_data, lines=lines, line_labels=line_labels, zlines=zoom_lines, zline_labels=zoom_line_labels),
        code="""
        var show_emission = false
        var show_absorption = false
        if (cb_obj.active.indexOf(0) >= 0) {  // index 0=Emission in active list
            show_emission = true
        }
        if (cb_obj.active.indexOf(1) >= 0) {  // index 1=Absorption in active list
            show_absorption = true
        }

        for(var i=0; i<lines.length; i++) {
            if(line_data.data['emission'][i]) {
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
    # lines_button_group.js_on_change('value', lines_callback)

    # edit visual inspection
    vi_commentinput = TextInput(value='-', title="VI comment :")
    vi_nameinput = TextInput(value=username, title="Your name :")
    viflags = ["Yes","No","Maybe","LowSNR","Bad"] # To put somewhere else in code
    #vi_flaginput = Select(title="VI flag :", value="-", options=viflags)
    vi_flaginput = RadioButtonGroup(labels=viflags)

    add_viflag_callback = CustomJS(args=dict(cds_targetinfo=cds_targetinfo,ifiberslider = ifiberslider, vi_flaginput=vi_flaginput, viflags=viflags, nspec=nspec), code="""
        cds_targetinfo.data['VI_ongoing_flag'][ifiberslider.value]=viflags[vi_flaginput.active]
   //     if(ifiberslider.value<nspec-1) {
   //        ifiberslider.value += 1 ;
   //     }
   //     vi_flaginput.active = -1 // Reset once scan is done
    """)
    add_vicomment_callback = CustomJS(args=dict(cds_targetinfo=cds_targetinfo,ifiberslider = ifiberslider, vi_commentinput=vi_commentinput), code="""
        cds_targetinfo.data['VI_ongoing_comment'][ifiberslider.value]=vi_commentinput.value
    """)
    change_viname_callback = CustomJS(args=dict(cds_targetinfo=cds_targetinfo,nspec = nspec, vi_nameinput=vi_nameinput), code="""
        for (var i=0; i<nspec; i++) {
            cds_targetinfo.data['VI_ongoing_scanner'][i]=vi_nameinput.value
        }
    """)
    vi_commentinput.js_on_change('value',add_vicomment_callback)
    vi_nameinput.js_on_change('value',change_viname_callback)
    vi_flaginput.js_on_click(add_viflag_callback)

    # save VI info to ASCII file
    # tested briefly safari chrome firefox
    # Warning text output very sensitve for # " \  ... (standard js formatting not ok)
    save_vi_button = Button(label="Download VI",button_type="default")
    save_vi_callback = CustomJS(args=dict(cds_targetinfo=cds_targetinfo, viflags=viflags), code="""
        function download(filename, text) {
            var element = document.createElement('a');
            element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(text));
            element.setAttribute('download', filename);
            element.style.display = 'none';
            document.body.appendChild(element);
            element.click();
            document.body.removeChild(element);
        }
        var thetext="# Prototype VI result file generated by prospect \\n";
        thetext+= "# TargetID Expid Fiber Spec_version Redrock_version Redrock_spectype Redrock_z VI_Scanner VI_flag VI_comment \\n";
        var toto = cds_targetinfo.data['VI_ongoing_flag'];
        var titi = cds_targetinfo.data['VI_ongoing_scanner'];
        var tutu = cds_targetinfo.data['VI_ongoing_comment'];
        var fa = cds_targetinfo.data['expid'];
        var fb = cds_targetinfo.data['fiber'];
        var fc = cds_targetinfo.data['targetid'];
        var fd = cds_targetinfo.data['spec_version'];
        var fe = cds_targetinfo.data['redrock_version'];
        var za = cds_targetinfo.data['spectype'] ;
        var zb = cds_targetinfo.data['z'] ;
        for (var i=0 ; i< toto.length; i++) {
            if (viflags.includes(toto[i])) thetext += (fc[i]+" "+fa[i]+" "+fb[i]+" "+fd[i]+" "+fe[i]+" "+za[i]+" "+zb[i].toFixed(3)+" "+titi[i]+" "+toto[i]+' "'+tutu[i]+'"'+" \\n");
        }
        download("vi_result.txt",thetext) ;
    """)
    save_vi_button.js_on_event('button_click', save_vi_callback)

    # Choose to show or not previous VI
    show_prev_vi_select = Select(title='Show previous VI', value='No', options=['Yes','No'])
    show_prev_vi_callback = CustomJS(args=dict(vi_info_div = vi_info_div, show_prev_vi_select=show_prev_vi_select, targetinfo = cds_targetinfo, ifiberslider = ifiberslider), code="""
        if (show_prev_vi_select.value == "Yes") {
            vi_info_div.text = targetinfo.data['vi_info'][ifiberslider.value];
        } else {
            vi_info_div.text = " ";
        }
    """)
    show_prev_vi_select.js_on_change('value',show_prev_vi_callback)

    
    #-----
    reset_plotrange_button = Button(label="Reset X-Y range",button_type="default")
    reset_plotrange_callback = CustomJS(args = dict(fig=fig, xmin=xmin, xmax=xmax, spectra=cds_spectra), code="""
        // x-range : use fixed x-range determined once for all
        fig.x_range.start = xmin
        fig.x_range.end = xmax
        
        // y-range : same function as in update_plot()
        function get_y_minmax(pmin, pmax, data) {
            var dx = data.slice().filter(Boolean)
            dx.sort()
            var imin = Math.floor(pmin * dx.length)
            var imax = Math.floor(pmax * dx.length)
            return [dx[imin], dx[imax]]
        }
        var ymin = 0.0
        var ymax = 0.0
        for (var i=0; i<spectra.length; i++) {
            var data = spectra[i].data
            tmp = get_y_minmax(0.01, 0.99, data['plotflux'])
            ymin = Math.min(ymin, tmp[0])
            ymax = Math.max(ymax, tmp[1])
        }
        if(ymin<0) {
            fig.y_range.start = ymin * 1.4
        } else {
            fig.y_range.start = ymin * 0.6
        }
        fig.y_range.end = ymax * 1.4

    """)
    reset_plotrange_button.js_on_event('button_click', reset_plotrange_callback)
    
    #-----
    update_plot = CustomJS(
        args = dict(
            spectra = cds_spectra,
            coaddcam_spec = cds_coaddcam_spec,
            model = cds_model,
            targetinfo = cds_targetinfo,
            target_info_div = target_info_div,
            vi_info_div = vi_info_div,
            show_prev_vi_select = show_prev_vi_select,
            ifiberslider = ifiberslider,
            smootherslider = smootherslider,
            zslider=zslider,
            dzslider=dzslider,
  #          lines_button_group = lines_button_group,
            fig = fig,
            imfig_source=imfig_source,
            imfig_urls=imfig_urls,
            vi_commentinput=vi_commentinput,
            vi_nameinput=vi_nameinput,
            vi_flaginput=vi_flaginput,
            viflags = viflags
            ),
        code = """
        var ifiber = ifiberslider.value
        var nsmooth = smootherslider.value
        target_info_div.text = targetinfo.data['target_info'][ifiber]
        if (show_prev_vi_select.value == "Yes") {
            vi_info_div.text = targetinfo.data['vi_info'][ifiber];
        } else {
            vi_info_div.text = " ";
        }
        // Added EA : ongoing VI
        if (cb_obj == ifiberslider) {
            vi_commentinput.value=targetinfo.data['VI_ongoing_comment'][ifiber] ;
            vi_nameinput.value=targetinfo.data['VI_ongoing_scanner'][ifiber] ;
         //   vi_flaginput.value=targetinfo.data['VI_ongoing_flag'][ifiber] ;
            vi_flaginput.active = viflags.indexOf(targetinfo.data['VI_ongoing_flag'][ifiber]) ; // -1 if nothing
        }

        if(targetinfo.data['z'] != undefined && cb_obj == ifiberslider) {
            var z = targetinfo.data['z'][ifiber]
            var z1 = Math.floor(z*100) / 100
            zslider.value = z1
            dzslider.value = (z - z1)
        }

        function get_y_minmax(pmin, pmax, data) {
            // copy before sorting to not impact original, and filter out NaN
            var dx = data.slice().filter(Boolean)
            dx.sort()
            var imin = Math.floor(pmin * dx.length)
            var imax = Math.floor(pmax * dx.length)
            return [dx[imin], dx[imax]]
        }

        // Smoothing kernel
        if (nsmooth > 0) {
            var kernel = [];
            for(var i=-2*nsmooth; i<=2*nsmooth; i++) {
                kernel.push(Math.exp(-(i**2)/(2*nsmooth)))
            }
            var kernel_offset = Math.floor(kernel.length/2)
        }

        function smooth_data(data_in, kernel, kernel_offset, quadrature=false) {
            // by default : out_j ~ (sum K_i in_i) / (sum K_i)
            // if quadrature is true (for noise) : out_j^2 ~ (sum K_i^2 in_i^2) / (sum K_i)^2
            var smoothed_data = data_in.slice()
            for (var j=0; j<data_in.length; j++) {
                smoothed_data[j] = 0.0
                var weight = 0.0
                // TODO: speed could be improved by moving `if` out of loop
                for (var k=0; k<kernel.length; k++) {
                    var m = j+k-kernel_offset
                    if((m >= 0) && (m < data_in.length)) {
                        var fx = data_in[m]
                        if(fx == fx) {
                            if (quadrature==true) {
                                smoothed_data[j] = smoothed_data[j] + (fx * kernel[k])**2
                            } else {
                                smoothed_data[j] = smoothed_data[j] + fx * kernel[k]
                            }
                            weight += kernel[k]
                        }
                    }
                }
                if (quadrature==true) {
                    smoothed_data[j] = Math.sqrt(smoothed_data[j]) / weight
                } else {
                    smoothed_data[j] = smoothed_data[j] / weight
                }
            }
            return smoothed_data
        }
        
        // Find nearest index in grid, left from point; use dichotomy method
        function index_dichotomy(point, grid) {
            if ( point < grid[0] ) return 0
            if ( point > grid[grid.length-1] ) return grid.length-2
            var i_left = 0
            var i_center = 0
            var i_right = grid.length-1
            while ( i_right - i_left != 1) {
                i_center = i_left + Math.floor((i_right-i_left)/2)
                if ( point >= grid[i_center] ) {
                    i_left = i_center
                } else {
                    i_right = i_center
                }
            }
            return i_left
        }

        // Basic linear interpolation at on point
        function interp_grid(xval, xarr, yarr) {
            var index = index_dichotomy(xval, xarr)
            var a = (yarr[index+1] - yarr[index])/(xarr[index+1] - xarr[index])
            var b = yarr[index]-a*xarr[index]
            var yval = a*xval+b
            return yval
        }

        // Coadd brz spectra. Similar to the python code mycoaddcam()
        function coadd_brz_cams(wave_in, flux_in, noise_in) {
            // each "_in" must have 3 entries (brz)
            // TODO handle case of no noise

            // Find b,r,z ordering in input arrays
            var wave_start = [wave_in[0][0], wave_in[1][0], wave_in[2][0]]
            var i_b = wave_start.indexOf(Math.min.apply(Math, wave_start))
            var i_z = wave_start.indexOf(Math.max.apply(Math, wave_start))
            var i_r = 1
            for (var i=0; i<3; i++) {
                if ( (i_b != i) && (i_z != i) ) i_r = i
            }

            var wave_out = []
            var flux_out = []
            var noise_out = []
            var margin = 20
            for (var i=0; i<wave_in[i_b].length; i++) { // b
                if (wave_in[i_b][i] < wave_in[i_b][wave_in[i_b].length-1] - margin) {
                    wave_out.push(wave_in[i_b][i])
                    flux_out.push(flux_in[i_b][i])
                    noise_out.push(noise_in[i_b][i])
                }
            }
            var the_lim = wave_out[wave_out.length-1]
            for (var i=0; i<wave_in[i_r].length; i++) { // r
                if ( (wave_in[i_r][i] < wave_in[i_r][wave_in[i_r].length-1] - margin) && (wave_in[i_r][i] > the_lim)) {
                    wave_out.push(wave_in[i_r][i])
                    flux_out.push(flux_in[i_r][i])
                    noise_out.push(noise_in[i_r][i])
                }
            }
            the_lim = wave_out[wave_out.length-1]
            for (var i=0; i<wave_in[i_z].length; i++) { // z
                if (wave_in[i_z][i] > the_lim) {
                    wave_out.push(wave_in[i_z][i])
                    flux_out.push(flux_in[i_z][i])
                    noise_out.push(noise_in[i_z][i])
                }
            }
            for (var i=0; i<wave_out.length; i++) { // combine in overlapping regions
                var b1 = -1
                var b2 = -1
                if ( (wave_out[i] > wave_in[i_r][0]) && (wave_out[i] < wave_in[i_b][wave_in[i_b].length-1]) ) { // br
                    b1 = 0
                    b2 = 1
                }
                if ( (wave_out[i] > wave_in[i_z][0]) && (wave_out[i] < wave_in[i_r][wave_in[i_r].length-1]) ) {  // rz
                    b1 = 1
                    b2 = 2
                }
                if (b1 != -1) {
                    var phi1 = interp_grid(wave_out[i], wave_in[b1], flux_in[b1])
                    var noise1 = interp_grid(wave_out[i], wave_in[b1], noise_in[b1])
                    var phi2 = interp_grid(wave_out[i], wave_in[b2], flux_in[b2])
                    var noise2 = interp_grid(wave_out[i], wave_in[b2], noise_in[b2])
                    if ( noise1 > 0 && noise2 > 0 ) {
                        var iv1 = 1/(noise1*noise1)
                        var iv2 = 1/(noise2*noise2)
                        var iv = iv1+iv2
                        noise_out[i] = 1/Math.sqrt(iv)
                        flux_out[i] = (iv1*phi1+iv2*phi2)/iv
                    }
                }
            }
            return [wave_out, flux_out, noise_out]
        }
        
        // Smooth plot and recalculate ymin/ymax
        var ymin = 0.0
        var ymax = 0.0
        for (var i=0; i<spectra.length; i++) {
            var data = spectra[i].data
            var origflux = data['origflux'+ifiber]
            if ('plotnoise' in data) {
                var orignoise = data['orignoise'+ifiber]
            }
            if (nsmooth == 0) {
                data['plotflux'] = origflux.slice()
                if ('plotnoise' in data) {
                    data['plotnoise'] = orignoise.slice()
                }
            } else {
                data['plotflux'] = smooth_data(origflux, kernel, kernel_offset)
                if ('plotnoise' in data) {
                    // Add noise in quadrature
                    data['plotnoise'] = smooth_data(orignoise, kernel, kernel_offset, quadrature=true)
                }
            }
            spectra[i].change.emit()

            tmp = get_y_minmax(0.01, 0.99, data['plotflux'])
            ymin = Math.min(ymin, tmp[0])
            ymax = Math.max(ymax, tmp[1])
        }

        // update camera-coadd
        // Here I choose to do coaddition on the smoothed spectra (should be ok?)
        var wave_in = []
        var flux_in = []
        var noise_in = []
        for (var i=0; i<3; i++) {
            var data = spectra[i].data
            wave_in.push(data['plotwave'].slice())
            flux_in.push(data['plotflux'].slice())
            if ('plotnoise' in data) {
                noise_in.push(data['plotnoise'].slice())
            } else {
                var dummy_noise = []
                for (var j=0; j<data['plotflux'].length; j++) dummy_noise.push(1)
                noise_in.push(dummy_noise)
            }
        }
        var coadd_infos = coadd_brz_cams(wave_in, flux_in, noise_in)
        coaddcam_spec.data['plotwave'] = coadd_infos[0].slice()
        coaddcam_spec.data['plotflux'] = coadd_infos[1].slice()
        coaddcam_spec.data['plotnoise'] = coadd_infos[2].slice()
        coaddcam_spec.change.emit()

        // update model
        if(model) {
            var origflux = model.data['origflux'+ifiber]
            if (nsmooth == 0) {
                model.data['plotflux'] = origflux.slice()
            } else {
                model.data['plotflux'] = smooth_data(origflux, kernel, kernel_offset)
            }
            model.change.emit()
        }

        // update y_range
        if(ymin<0) {
            fig.y_range.start = ymin * 1.4
        } else {
            fig.y_range.start = ymin * 0.6
        }
        fig.y_range.end = ymax * 1.4
        //
        // update target image
        //
        imfig_source.data.url[0] = imfig_urls[ifiber][0]
        imfig_source.data.txt[0] = imfig_urls[ifiber][2]
        imfig_source.change.emit()
    """)
    smootherslider.js_on_change('value', update_plot)
    ifiberslider.js_on_change('value', update_plot)

    #-----
    #- Add navigation buttons
    navigation_button_width = 30
    prev_button = Button(label="<", width=navigation_button_width)
    next_button = Button(label=">", width=navigation_button_width)
    
    prev_callback = CustomJS(
        args=dict(ifiberslider=ifiberslider),
        code="""
        if(ifiberslider.value>0) {
            ifiberslider.value--
        }
        """)
    next_callback = CustomJS(
        args=dict(ifiberslider=ifiberslider, nspec=nspec),
        code="""
        if(ifiberslider.value<nspec-1) {
            ifiberslider.value++
        }
        """)

    prev_button.js_on_event('button_click', prev_callback)
    next_button.js_on_event('button_click', next_callback)
    
    #-----
    slider_width = plot_width - 2*navigation_button_width
    navigator = bk.Row(
        widgetbox(prev_button, width=navigation_button_width+15),
        widgetbox(vi_flaginput, width=60*len(viflags)),
        widgetbox(next_button, width=navigation_button_width+20),
        widgetbox(ifiberslider, width=plot_width-(60*len(viflags)+2*navigation_button_width+40)))
    the_bokehsetup = bk.Column(
            bk.Row(fig, bk.Column(imfig, zoomfig)),
            bk.Row(
                widgetbox(target_info_div, width=plot_width - 120),
                widgetbox(reset_plotrange_button, width = 100)
            ),
            navigator,
            bk.Row(
                widgetbox(smootherslider, width=plot_width//2),
                widgetbox(display_options_group,width=120),
            ),
            bk.Row(
                widgetbox(waveframe_buttons, width=120),
                widgetbox(zslider, width=plot_width//2 - 120),
                widgetbox(dzslider, width=plot_width//2 - 120),
                widgetbox(z_display, width=120),
                widgetbox(zreset_button, width=100)
            ),
            bk.Row(
                widgetbox(lines_button_group),
                widgetbox(coaddcam_buttons, width=120)
            ),
            bk.Row(
                widgetbox(vi_commentinput,width=plot_width-500),
                widgetbox(vi_nameinput,width=120),
                widgetbox(save_vi_button,width=100,sizing_mode="scale_height")
            ),
#            widgetbox(save_vi_button,width=100)
            ## Don't want this in principle :
#            bk.Row(
#                widgetbox(show_prev_vi_select,width=100),
#                widgetbox(vi_info_div, width=plot_width-130)
#                )
        )
    if notebook:
        bk.show(the_bokehsetup)
    else:
        bk.save(the_bokehsetup)
    

#-------------------------------------------------------------------------
_line_list = [
    #
    # This is the set of emission lines from the spZline files.
    # See $IDLSPEC2D_DIR/etc/emlines.par
    # Wavelengths are in air for lambda > 2000, vacuum for lambda < 2000.
    # TODO: convert to vacuum wavelengths
    #
    {"name" : "Lyα",      "longname" : "Lyman α",        "lambda" : 1215.67,  "emission": True },
    {"name" : "Lyβ",      "longname" : "Lyman β",        "lambda" : 1025.18,  "emission": True },
    {"name" : "N V",      "longname" : "N V 1240",       "lambda" : 1240.81,  "emission": True },
    {"name" : "C IV",     "longname" : "C IV 1549",      "lambda" : 1549.48,  "emission": True },
    {"name" : "He II",    "longname" : "He II 1640",     "lambda" : 1640.42,  "emission": True },
    {"name" : "C III]",   "longname" : "C III] 1908",    "lambda" : 1908.734, "emission": True },
    {"name" : "Mg II",    "longname" : "Mg II 2799",     "lambda" : 2799.49,  "emission": True },
    {"name" : "[O II]",   "longname" : "[O II] 3725",    "lambda" : 3726.032, "emission": True },
    {"name" : "[O II]",   "longname" : "[O II] 3727",    "lambda" : 3728.815, "emission": True },
    {"name" : "[Ne III]", "longname" : "[Ne III] 3868",  "lambda" : 3868.76,  "emission": True },
    {"name" : "Hζ",       "longname" : "Balmer ζ",       "lambda" : 3889.049, "emission": True },
    {"name" : "[Ne III]", "longname" : "[Ne III] 3970",  "lambda" : 3970.00,  "emission": True },
    {"name" : "Hε",       "longname" : "Balmer ε",       "lambda" : 3970.072, "emission": True },
    {"name" : "Hδ",       "longname" : "Balmer δ",       "lambda" : 4101.734, "emission": True },
    {"name" : "Hγ",       "longname" : "Balmer γ",       "lambda" : 4340.464, "emission": True },
    {"name" : "[O III]",  "longname" : "[O III] 4363",   "lambda" : 4363.209, "emission": True },
    {"name" : "He II",    "longname" : "He II 4685",     "lambda" : 4685.68,  "emission": True },
    {"name" : "Hβ",       "longname" : "Balmer β",       "lambda" : 4861.325, "emission": True },
    {"name" : "[O III]",  "longname" : "[O III] 4959",   "lambda" : 4958.911, "emission": True },
    {"name" : "[O III]",  "longname" : "[O III] 5007",   "lambda" : 5006.843, "emission": True },
    {"name" : "He II",    "longname" : "He II 5411",     "lambda" : 5411.52,  "emission": True },
    {"name" : "[O I]",    "longname" : "[O I] 5577",     "lambda" : 5577.339, "emission": True },
    {"name" : "[N II]",   "longname" : "[N II] 5755",    "lambda" : 5754.59,  "emission": True },
    {"name" : "He I",     "longname" : "He I 5876",      "lambda" : 5875.68,  "emission": True },
    {"name" : "[O I]",    "longname" : "[O I] 6300",     "lambda" : 6300.304, "emission": True },
    {"name" : "[S III]",  "longname" : "[S III] 6312",   "lambda" : 6312.06,  "emission": True },
    {"name" : "[O I]",    "longname" : "[O I] 6363",     "lambda" : 6363.776, "emission": True },
    {"name" : "[N II]",   "longname" : "[N II] 6548",    "lambda" : 6548.05,  "emission": True },
    {"name" : "Hα",       "longname" : "Balmer α",       "lambda" : 6562.801, "emission": True },
    {"name" : "[N II]",   "longname" : "[N II] 6583",    "lambda" : 6583.45,  "emission": True },
    {"name" : "[S II]",   "longname" : "[S II] 6716",    "lambda" : 6716.44,  "emission": True },
    {"name" : "[S II]",   "longname" : "[S II] 6730",    "lambda" : 6730.82,  "emission": True },
    {"name" : "[Ar III]", "longname" : "[Ar III] 7135",  "lambda" : 7135.790, "emission": True },
    #
    # Absorption lines
    #
    {"name" : "Hζ",   "longname" : "Balmer ζ",         "lambda" : 3889.049, "emission": False },
    {"name" : "K",    "longname" : "K (Ca II 3933)",   "lambda" : 3933.7,   "emission": False },
    {"name" : "H",    "longname" : "H (Ca II 3968)",   "lambda" : 3968.5,   "emission": False },
    {"name" : "Hε",   "longname" : "Balmer ε",         "lambda" : 3970.072, "emission": False },
    {"name" : "Hδ",   "longname" : "Balmer δ",         "lambda" : 4101.734, "emission": False },
    {"name" : "G",    "longname" : "G (Ca I 4307)",    "lambda" : 4307.74,  "emission": False },
    {"name" : "Hγ",   "longname" : "Balmer γ",         "lambda" : 4340.464, "emission": False },
    {"name" : "Hβ",   "longname" : "Balmer β",         "lambda" : 4861.325, "emission": False },
    {"name" : "Mg I", "longname" : "Mg I 5175",        "lambda" : 5175.0,   "emission": False },
    {"name" : "D2",   "longname" : "D2 (Na I 5889)",   "lambda" : 5889.95,  "emission": False },
    # {"name" : "D",    "longname" : "D (Na I doublet)", "lambda": 5892.9,   "emission": False },
    {"name" : "D1",   "longname" : "D1 (Na I 5895)",   "lambda" : 5895.92,  "emission": False },
    {"name" : "Hα",   "longname" : "Balmer α",         "lambda" : 6562.801, "emission": False },
  ]

def _airtovac(w):
    """Convert air wavelengths to vacuum wavelengths. Don't convert less than 2000 Å.

    Parameters
    ----------
    w : :class:`float`
        Wavelength [Å] of the line in air.

    Returns
    -------
    :class:`float`
        Wavelength [Å] of the line in vacuum.
    """
    if w < 2000.0:
        return w;
    vac = w
    for iter in range(2):
        sigma2 = (1.0e4/vac)*(1.0e4/vac)
        fact = 1.0 + 5.792105e-2/(238.0185 - sigma2) + 1.67917e-3/(57.362 - sigma2)
        vac = w*fact
    return vac

def add_lines(fig, z=0 , emission=True, fig_height=None, label_offsets=[100, 5]):
    """
    label_offsets = [offset_absorption_lines, offset_emission_lines] : offsets in y-position 
                    for line labels wrt top (resp. bottom) of the figure
    """
    
    if fig_height is None : fig_height = fig.plot_height

    line_data = dict()
    line_data['restwave'] = np.array([_airtovac(row['lambda']) for row in _line_list])
    line_data['plotwave'] = line_data['restwave'] * (1+z)
    line_data['name'] = [row['name'] for row in _line_list]
    line_data['longname'] = [row['name'] for row in _line_list]
    line_data['plotname'] = [row['name'] for row in _line_list]
    line_data['emission'] = [row['emission'] for row in _line_list]

    y = list()
    for i in range(len(line_data['restwave'])):
        if i == 0:
            if _line_list[i]['emission']:
                y.append(fig_height - label_offsets[0])
            else:
                y.append(label_offsets[1])
        else:
            if (line_data['restwave'][i] < line_data['restwave'][i-1]+label_offsets[0]) and \
               (line_data['emission'][i] == line_data['emission'][i-1]):
                if line_data['emission'][i]:
                    y.append(y[-1] - 15)
                else:
                    y.append(y[-1] + 15)
            else:
                if line_data['emission'][i]:
                    y.append(fig_height-label_offsets[0])
                else:
                    y.append(label_offsets[1])

    line_data['y'] = y

    #- Add vertical spans to figure
    lines = list()
    labels = list()
    for w, y, name, emission in zip(
            line_data['plotwave'],
            line_data['y'],
            line_data['plotname'],
            line_data['emission']
            ):
        if emission:
            color = 'magenta'
        else:
            color = 'green'

        s = Span(location=w, dimension='height', line_color=color,
                line_alpha=1.0, visible=False)

        fig.add_layout(s)
        lines.append(s)

        lb = Label(x=w, y=y, x_units='data', y_units='screen',
                    text=name, text_color='gray', text_font_size="8pt",
                    x_offset=2, y_offset=0, visible=False)
        fig.add_layout(lb)
        labels.append(lb)

    line_data = bk.ColumnDataSource(line_data)
    return line_data, lines, labels


if __name__ == '__main__':
    # framefiles = [
    #     'data/cframe-b0-00000020.fits',
    #     'data/cframe-r0-00000020.fits',
    #     'data/cframe-z0-00000020.fits',
    # ]
    #
    # frames = list()
    # for filename in framefiles:
    #     fr = desispec.io.read_frame(filename)
    #     fr = fr[0:50]  #- Trim for faster debugging
    #     frames.append(fr)
    #
    # plotspectra(frames)

    # Outdated :

    parser = argparse.ArgumentParser(description='Create html pages for the spectral viewer')
    parser.add_argument('healpixel', help='Healpixel (nside64) to process', type=str)
    parser.add_argument('--basedir', help='Path to spectra reltive to DESI_ROOT', type=str, default="datachallenge/reference_runs/18.6/spectro/redux/mini/spectra-64")
    args = parser.parse_args()
    basedir = os.environ['DESI_ROOT']+"/"+args.basedir+"/"+args.healpixel[0:2]+"/"+args.healpixel+"/"

    specfile = basedir+'spectra-64-'+args.healpixel+'.fits'
    zbfile = specfile.replace('spectra-64-', 'zbest-64-')

    #- Original remapping of individual spectra to zbest
    # spectra = desispec.io.read_spectra(specfile)
    # zbest_raw = Table.read(zbfile, 'ZBEST')

    # # EA : all is best is zbest matches spectra row-by-row.
    # zbest=Table(dtype=zbest_raw.dtype)
    # for i in range(spectra.num_spectra()) :
    #     ww, = np.where((zbest_raw['TARGETID'] == spectra.fibermap['TARGETID'][i]))
    #     if len(ww)!=1 : print("!! Issue with zbest table !!")
    #     zbest.add_row(zbest_raw[ww[0]])

    #- Coadd on the fly
    individual_spectra = desispec.io.read_spectra(specfile)
    spectra = coadd_targets(individual_spectra)
    zbest = Table.read(zbfile, 'ZBEST')

    mwave, mflux = create_model(spectra, zbest)

    ## VI "catalog" - location to define later
    vifile = os.environ['HOME']+"/prospect/vilist_prototype.fits"
    vidata = utils_specviewer.match_vi_targets(vifile, spectra.fibermap['TARGETID'])

    plotspectra(spectra, zcatalog=zbest, vidata=vidata, model=(mwave, mflux), title=os.path.basename(specfile))

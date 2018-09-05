import os, sys

import numpy as np
from astropy.table import Table
import bokeh.plotting as bk

from bokeh.models import ColumnDataSource, CDSView, IndexFilter
from bokeh.models import CustomJS, LabelSet, Label, Span
from bokeh.models.widgets import (
    Slider, Button, Div, CheckboxButtonGroup, RadioButtonGroup)
from bokeh.layouts import widgetbox
import bokeh.events
# from bokeh.layouts import row, column

import desispec.io
from desitarget.targetmask import desi_mask
import desispec.spectra
import desispec.frame

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


def plotspectra(spectra, zcatalog=None, model=None, notebook=False, title=None):
    '''
    TODO: document
    '''

    if notebook:
        bk.output_notebook()

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

        cdsdata['plotflux'] = cdsdata['origflux0']

        cds_spectra.append(
            bk.ColumnDataSource(cdsdata, name=band)
            )

    #- Gather models into ColumnDataSource objects, row matched to spectra
    #- TODO: allow more than one zcatalog entry with different ZNUM per targetid
    mwave, mflux = model
    model_obswave = mwave.copy()
    model_restwave = mwave.copy()
    cds_model_data = dict(
        origwave = mwave.copy(),
        plotwave = mwave.copy(),
        plotflux = np.zeros(len(mwave)),
    )
    targetids = spectra.target_ids()
    for i in range(nspec):
        key = 'origflux'+str(i)
        j = np.where(zcatalog['TARGETID'] == targetids[i])[0][0]
        cds_model_data[key] = mflux[j]

    cds_model_data['plotflux'] = cds_model_data['origflux0']

    cds_model = bk.ColumnDataSource(cds_model_data)

    #- TODO: add redshifts here?
    fmdict = dict()
    for colname in ['TARGETID', 'DESI_TARGET']:
        fmdict[colname] = spectra.fibermap[colname]

    target_names = list()
    for dt in spectra.fibermap['DESI_TARGET']:
        names = ' '.join(desi_mask.names(dt))
        target_names.append(names)

    fmdict['TARGET_NAMES'] = target_names
    cds_fibermap = bk.ColumnDataSource(fmdict, name='fibermap')

    plot_width=800
    plot_height=400
    # tools = 'pan,box_zoom,wheel_zoom,undo,redo,reset,save'
    tools = 'pan,box_zoom,wheel_zoom,reset,save'
    fig = bk.figure(height=plot_height, width=plot_width, title=title,
        tools=tools, toolbar_location='above', y_range=(-10, 20))
    fig.toolbar.active_drag = fig.tools[1]    #- box zoom
    fig.toolbar.active_scroll = fig.tools[2]  #- wheel zoom
    fig.xaxis.axis_label = 'Wavelength [Å]'
    fig.yaxis.axis_label = 'Flux'
    fig.xaxis.axis_label_text_font_style = 'normal'
    fig.yaxis.axis_label_text_font_style = 'normal'
    colors = dict(b='#1f77b4', r='#d62728', z='maroon')
    for spec in cds_spectra:
        fig.line('plotwave', 'plotflux', source=spec, line_color=colors[spec.name])

    fig.line('plotwave', 'plotflux', source=cds_model, line_color='black')

    zoomfig = bk.figure(height=plot_height//2, width=plot_height//2,
        y_range=fig.y_range, x_range=(5000,5100),
        # output_backend="webgl",
        toolbar_location=None, tools=[])

    for spec in cds_spectra:
        zoomfig.line('plotwave', 'plotflux', source=spec,
            line_color=colors[spec.name], line_width=1, line_alpha=1.0)

    zoomfig.line('plotwave', 'plotflux', source=cds_model, line_color='black')

    #- Callback to update zoom window x-range
    zoom_callback = CustomJS(
        args=dict(zoomfig=zoomfig),
        code="""
            zoomfig.x_range.start = cb_obj.x - 100;
            zoomfig.x_range.end = cb_obj.x + 100;
        """)

    fig.js_on_event(bokeh.events.MouseMove, zoom_callback)

    #-----
    #- Emission and absorption lines
    line_data, lines, line_labels = add_lines(fig, z=0)

    #-----
    #- Add widgets for controling plots
    zslider = Slider(start=0.0, end=4.0, value=0.0, step=0.01, title='Redshift')
    dzslider = Slider(start=0.0, end=0.01, value=0.0, step=0.0001, title='+ Delta redshift')
    dzslider.format = "0[.]0000"

    #- Observer vs. Rest frame wavelengths
    waveframe_buttons = RadioButtonGroup(
        labels=["Obs", "Rest"], active=0)

    zslider_callback  = CustomJS(
        args=dict(
            spectra=cds_spectra,
            model=cds_model,
            zslider=zslider,
            dzslider=dzslider,
            waveframe_buttons=waveframe_buttons,
            line_data=line_data, lines=lines, line_labels=line_labels,
            fig=fig,
            ),
        #- TODO: reorder to reduce duplicated code
        code="""
        var z = zslider.value + dzslider.value
        var line_restwave = line_data.data['restwave']
        // Observer Frame
        if(waveframe_buttons.active == 0) {
            var x = 0.0
            for(var i=0; i<line_restwave.length; i++) {
                x = line_restwave[i] * (1+z)
                lines[i].location = x
                line_labels[i].x = x
            }
            for(var i=0; i<spectra.length; i++) {
                var data = spectra[i].data
                var origwave = data['origwave']
                var plotwave = data['plotwave']
                for (var j=0; j<plotwave.length; j++) {
                    plotwave[j] = origwave[j]
                }
                spectra[i].change.emit()
            }

            // Update model wavelength array
            var origwave = model.data['origwave']
            var plotwave = model.data['plotwave']
            for(var i=0; i<plotwave.length; i++) {
                plotwave[i] = origwave[i]
            }
            model.change.emit()

        // Rest Frame
        } else {
            for(i=0; i<line_restwave.length; i++) {
                lines[i].location = line_restwave[i]
                line_labels[i].x = line_restwave[i]
            }
            for (var i=0; i<spectra.length; i++) {
                var data = spectra[i].data
                var origwave = data['origwave']
                var plotwave = data['plotwave']
                for (var j=0; j<plotwave.length; j++) {
                    plotwave[j] = origwave[j] / (1+z)
                }
                spectra[i].change.emit()
            }

            // Update model wavelength array
            var origwave = model.data['origwave']
            var plotwave = model.data['plotwave']
            for(var i=0; i<plotwave.length; i++) {
                plotwave[i] = origwave[i] / (1+z)
            }
            model.change.emit()
        }
        """)

    zslider.js_on_change('value', zslider_callback)
    dzslider.js_on_change('value', zslider_callback)
    waveframe_buttons.js_on_click(zslider_callback)

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

    ifiberslider = Slider(start=0, end=nspec-1, value=0, step=1)
    if frame_input:
        ifiberslider.title = 'Fiber'
    else:
        ifiberslider.title = 'Target'

    smootherslider = Slider(start=1, end=31, value=1, step=1, title='Smooth')
    target_info = Div(text=target_names[0])

    #-----
    #- Toggle lines
    lines_button_group = CheckboxButtonGroup(
            labels=["Emission", "Absorption"], active=[])

    lines_callback = CustomJS(
        args = dict(line_data=line_data, lines=lines, line_labels=line_labels),
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
            } else {
                lines[i].visible = show_absorption
                line_labels[i].visible = show_absorption
            }
        }
        """
    )
    lines_button_group.js_on_click(lines_callback)
    # lines_button_group.js_on_change('value', lines_callback)

    #-----
    update_plot = CustomJS(
        args = dict(
            spectra = cds_spectra,
            model = cds_model,
            fibermap = cds_fibermap,
            ifiberslider = ifiberslider,
            smootherslider = smootherslider,
            lines_button_group = lines_button_group,
            target_info = target_info,
            fig = fig,
            ),
        #- TODO: add smoother function to reduce duplicated code
        code = """
        var ifiber = ifiberslider.value
        var nsmooth = smootherslider.value
        var target_names = fibermap.data['TARGET_NAMES']
        target_info.text = target_names[ifiber]

        function get_y_minmax(pmin, pmax, data) {
            // copy before sorting to not impact original, and filter out NaN
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
            var plotflux = data['plotflux']
            var origflux = data['origflux'+ifiber]
            for (var j=0; j<plotflux.length; j++) {
                plotflux[j] = 0.0
                var n = 0
                for (var k=Math.max(0, j-nsmooth); k<Math.min(plotflux.length, j+nsmooth); k++) {
                    fx = origflux[k]
                    if(fx == fx) {
                        plotflux[j] = plotflux[j] + fx
                        n++
                    }
                }
                plotflux[j] = plotflux[j] / n
            }
            spectra[i].change.emit()

            tmp = get_y_minmax(0.01, 0.99, plotflux)
            ymin = Math.min(ymin, tmp[0])
            ymax = Math.max(ymax, tmp[1])
        }

        // update model
        var plotflux = model.data['plotflux']
        var origflux = model.data['origflux'+ifiber]
        for (var j=0; j<plotflux.length; j++) {
            plotflux[j] = 0.0
            var n = 0
            for (var k=Math.max(0, j-nsmooth); k<Math.min(plotflux.length, j+nsmooth); k++) {
                fx = origflux[k]
                if(fx == fx) {
                    plotflux[j] = plotflux[j] + fx
                    n++
                }
            }
            plotflux[j] = plotflux[j] / n
        }
        model.change.emit()

        // update y_range
        if(ymin<0) {
            fig.y_range.start = ymin * 1.4
        } else {
            fig.y_range.start = ymin * 0.6
        }
        fig.y_range.end = ymax * 1.4
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
        if(ifiberslider.value<nspec+1) {
            ifiberslider.value++
        }
        """)

    prev_button.js_on_event('button_click', prev_callback)
    next_button.js_on_event('button_click', next_callback)

    #-----
    slider_width = plot_width - 2*navigation_button_width
    navigator = bk.Row(
        widgetbox(prev_button, width=navigation_button_width),
        widgetbox(next_button, width=navigation_button_width+20),
        widgetbox(ifiberslider, width=slider_width-20))
    bk.show(bk.Column(
        bk.Row(fig, zoomfig),
        widgetbox(target_info),
        navigator,
        widgetbox(smootherslider, width=plot_width//2),
        bk.Row(
            widgetbox(waveframe_buttons, width=120),
            widgetbox(zslider, width=plot_width//2 - 60),
            widgetbox(dzslider, width=plot_width//2 - 60),
            ),
        widgetbox(lines_button_group),
        ))
 
    #--- DEBUG ---
    # import IPython
    # IPython.embed()
    #--- DEBUG ---

#-------------------------------------------------------------------------
_line_list = [
    #
    # This is the set of emission lines from the spZline files.
    # See $IDLSPEC2D_DIR/etc/emlines.par
    # Wavelengths are in air for lambda > 2000, vacuum for lambda < 2000.
    # TODO: convert to vacuum wavelengths
    #
    {"name" : "Lyα",      "longname" : "Lyman α",        "lambda" : 1215.67,  "emission": True },
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

def add_lines(fig, z, emission=True, fig_height=350):
    line_data = dict()
    line_data['restwave'] = np.array([row['lambda'] for row in _line_list])
    line_data['plotwave'] = line_data['restwave'] * (1+z)
    line_data['name'] = [row['name'] for row in _line_list]
    line_data['longname'] = [row['name'] for row in _line_list]
    line_data['plotname'] = [row['name'] for row in _line_list]
    line_data['emission'] = [row['emission'] for row in _line_list]

    y = list()
    for i in range(len(line_data['restwave'])):
        if i == 0:
            if _line_list[i]['emission']:
                y.append(fig_height - 100)
            else:
                y.append(5)
        else:
            if (line_data['restwave'][i] < line_data['restwave'][i-1]+100) and \
               (line_data['emission'][i] == line_data['emission'][i-1]):
                if line_data['emission'][i]:
                    y.append(y[-1] - 15)
                else:
                    y.append(y[-1] + 15)
            else:
                if line_data['emission'][i]:
                    y.append(fig_height-100)
                else:
                    y.append(5)

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

    specfile = 'data/spectra-64-5261.fits'
    zbfile = specfile.replace('spectra-64-', 'zbest-64-')
    spectra = desispec.io.read_spectra(specfile)
    zbest = Table.read(zbfile, 'ZBEST')
    mwave, mflux = create_model(spectra, zbest)
    plotspectra(spectra, zcatalog=zbest, model=(mwave, mflux), title=os.path.basename(specfile))

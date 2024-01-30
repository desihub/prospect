# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
=====================
prospect.viewer.plots
=====================

Class containing bokeh plots needed for the viewer

"""

import numpy as np
_specutils_imported = True
try:
    from specutils import Spectrum1D, SpectrumList
except ImportError:
    _specutils_imported = False

import bokeh.plotting as bk
from bokeh.models import CustomJS, ColumnDataSource, BoxAnnotation, Legend, Span, Label
import bokeh.layouts as bl
import bokeh.events


def _cross_hair_points(x_cross, y_cross):
    '''
    Points to make a cross-hair centered on x_cross,y_cross
    Return (xs,ys) to be given to bokeh.plotting.figure.multi_line()
    '''
    xs = [[x_cross-15, x_cross-5], [x_cross+15, x_cross+5],
          [x_cross, x_cross],[x_cross, x_cross]]
    ys = [[y_cross, y_cross],[y_cross, y_cross],
          [y_cross-15, y_cross-5],[y_cross+5, y_cross+15]]
    return (xs, ys)


def _viewer_urls(spectra, pixscale=0.262, zoom=13, layer='ls-dr9'):
    """Return legacysurvey.org viewer URLs + other infos for the imaging cutout, for all spectra.

    Notes:
        `layer` does not apply to the JPEG cutout service.
        `pixscale` is in arcsec/pixel (default 0.262 is the native DECam scale)
    """

    #- Template for jpeg-cutout url
    u = "https://www.legacysurvey.org/viewer/jpeg-cutout?ra={0:f}&dec={1:f}&pixscale={2:f}"
    #- Template for on-click link to full viewer
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

    #- Compute PM corrections
    default_ref_epoch = 2015.5  # PM correction is set to zero at that epoch
    if hasattr(spectra, 'fibermap'):
        pmcor_ra = (default_ref_epoch-spectra.fibermap['REF_EPOCH'])*spectra.fibermap['PMRA']/1e3  # PM in mas/yr
        pmcor_dec = (default_ref_epoch-spectra.fibermap['REF_EPOCH'])*spectra.fibermap['PMDEC']/1e3
        # avoid adding a second cross-hair when PM correction is smaller than 1 arcsec (ie. in most cases):
        mask = ((np.abs(pmcor_ra)<1) & (np.abs(pmcor_dec)<1)) | (spectra.fibermap['REF_EPOCH']==0)
        pmcor_ra[mask] = 0
        pmcor_dec[mask] = 0
    else:
        pmcor_ra = spectra.meta['plugmap']['RA'] * 0.0
        pmcor_dec = spectra.meta['plugmap']['DEC'] * 0.0
    # convert PM correction to pixels in the jpeg cutout
    npix_cutout = 256
    x_pm = npix_cutout//2 - pmcor_ra/pixscale   # minus sign: RA decreases with x-coord in image
    y_pm = npix_cutout//2 + pmcor_dec/pixscale
    # return points to make a second cross-hair if PM correction is large
    list_crosshairs = [ _cross_hair_points(x_pm[i], y_pm[i]) for i in range(len(ra)) ]

    return [(u.format(ra[i], dec[i], pixscale),
             v.format(ra[i], dec[i], zoom, layer),
             'RA, Dec = {0:.4f}, {1:+.4f}'.format(ra[i], dec[i]),
             list_crosshairs[i][0], list_crosshairs[i][1])
            for i in range(len(ra))]


class ViewerPlots(object):
    """
    Encapsulates Bokeh plot-like objects that are part of prospect's GUI.
    """

    def __init__(self):
        # "Hardcoded" plotting parameters here:
        self.legend_outside_plot = True
        if (self.legend_outside_plot):
            self.xmargin_left = 100.
            self.xmargin_right = 100.
        else:
            self.xmargin_left = 200.
            self.xmargin_right = 400.
        self.plot_width=800
        self.plot_height=400
        self.colors = dict(b='#1f77b4', r='#d62728', z='maroon', coadd='#d62728', brz='#d62728')
        self.noise_colors = dict(b='greenyellow', r='green', z='forestgreen', coadd='green', brz='green')
        ## overlap wavelengths are hardcoded, from 1907.10688 (Table 1)
        self.overlap_waves = [ [5660, 5930], [7470, 7720] ]
        self.alpha_overlapband = 0.03

        self.fig = None
        self.zoomfig = None
        self.zoom_callback = None
        self.imfig = bl.Spacer(width=self.plot_height//2, height=self.plot_height//2)
        self.imfig_source = self.imfig_urls = self.crosshair_source = None


    def create_mainfig(self, spectra, title, viewer_cds, survey, with_noise=True, with_coaddcam=True):
        #-----
        #- Main figure
        #- Determine initial ymin, ymax, xmin, xmax
        self.ymin = self.ymax = self.xmax = 0
        self.xmin = 100000.
        if survey == 'SDSS':
            bands = ['coadd']
            try:
                self.ymin = np.nanmin(spectra.flux.value[0, :])
                self.ymax = np.nanmax(spectra.flux.value[0, :])
            except IndexError:
                # Catch case where, e.g., flux.shape == (3000, ).
                self.ymin = np.nanmin(spectra.flux.value)
                self.ymax = np.nanmax(spectra.flux.value)
            self.xmin = np.min(spectra.spectral_axis.value)
            self.xmax = np.max(spectra.spectral_axis.value)
        else:
            bands = spectra.bands
            for i, band in enumerate(bands):
                if _specutils_imported and isinstance(spectra, SpectrumList):
                    sp_flux = spectra[i].flux.value[0]
                    sp_wave = spectra[i].spectral_axis.value
                else:
                    sp_flux = spectra.flux[band][0]
                    sp_wave = spectra.wave[band]
                if np.isfinite(sp_flux).any():
                    self.ymin = min(self.ymin, np.nanmin(sp_flux))
                    self.ymax = max(self.ymax, np.nanmax(sp_flux))
                self.xmin = min(self.xmin, np.min(sp_wave))
                self.xmax = max(self.xmax, np.max(sp_wave))
        self.xmin -= self.xmargin_left
        self.xmax += self.xmargin_right

        tools = 'pan,box_zoom,ywheel_zoom,xwheel_zoom,save'
        tooltips_fig = [("wave","$x"),("flux","$y")]
        self.fig = bk.figure(height=self.plot_height, width=self.plot_width, title=title,
            tools=tools, toolbar_location='above', tooltips=tooltips_fig,
            y_range=(self.ymin, self.ymax), x_range=(self.xmin, self.xmax))
        self.fig.sizing_mode = 'stretch_width'
        self.fig.toolbar.active_drag = self.fig.tools[0]    #- pan zoom (previously box)
        self.fig.toolbar.active_scroll = self.fig.tools[2]  #- wheel zoom
        self.fig.xaxis.axis_label = 'Wavelength [Å]'
        self.fig.yaxis.axis_label = 'Flux [10⁻¹⁷ erg cm⁻² s⁻¹ Å⁻¹]'
        self.fig.xaxis.axis_label_text_font_style = 'normal'
        self.fig.yaxis.axis_label_text_font_style = 'normal'
        self.alpha_discrete = 0.2 # alpha for "almost-hidden" curves (single-arm spectra and noise by default)
        if not with_coaddcam : self.alpha_discrete = 1

        #- Highlight overlap regions between arms
        self.overlap_bands = []
        if bands == ['brz'] or set(bands) == set(['b','r','z']) :
            for i in range(len(self.overlap_waves)) :
                fill_alpha = self.alpha_overlapband # if with_coaddcam else 0
                self.overlap_bands.append( BoxAnnotation(left=self.overlap_waves[i][0], right=self.overlap_waves[i][1], fill_color='blue', fill_alpha=fill_alpha, line_alpha=0) )
                self.fig.add_layout(self.overlap_bands[-1])

        self.data_lines = list()
        for spec in viewer_cds.cds_spectra:
            lx = self.fig.line('plotwave', 'plotflux', source=spec,
                            line_color=self.colors[spec.name], line_alpha=self.alpha_discrete)
            self.data_lines.append(lx)
        if with_coaddcam :
            lx = self.fig.line('plotwave', 'plotflux', source=viewer_cds.cds_coaddcam_spec,
                            line_color=self.colors['coadd'], line_alpha=1)
            self.data_lines.append(lx)

        self.noise_lines = list()
        if with_noise :
            for spec in viewer_cds.cds_spectra :
                lx = self.fig.line('plotwave', 'plotnoise', source=spec,
                                line_color=self.noise_colors[spec.name], line_alpha=self.alpha_discrete)
                self.noise_lines.append(lx)
            if with_coaddcam :
                lx = self.fig.line('plotwave', 'plotnoise', source=viewer_cds.cds_coaddcam_spec,
                                line_color=self.noise_colors['coadd'], line_alpha=1)
                self.noise_lines.append(lx)

        self.model_lines = list()
        if viewer_cds.cds_model is not None:
            lx = self.fig.line('plotwave', 'plotflux', source=viewer_cds.cds_model, line_color='black')
            self.model_lines.append(lx)

        self.othermodel_lines = list()
        if viewer_cds.cds_othermodel is not None :
            lx = self.fig.line('plotwave', 'plotflux', source=viewer_cds.cds_othermodel, line_color='black', line_dash='dashed')
            self.othermodel_lines.append(lx)

        legend_items = [("data", self.data_lines[-1::-1])] #- reversed to get blue as lengend entry
        if viewer_cds.cds_model is not None :
            legend_items.append(("pipeline fit", self.model_lines))
        if viewer_cds.cds_othermodel is not None :
            legend_items.append(("other model", self.othermodel_lines))
        if with_noise :
            legend_items.append(("noise", self.noise_lines[-1::-1])) # same as for data_lines
        if (self.legend_outside_plot):
            legend = Legend(items=legend_items,
                            border_line_alpha=0, label_text_font_size='11px',
                            glyph_width=20, margin=0, padding=0, spacing=20,
                            orientation='horizontal')
            self.fig.add_layout(legend, 'below')
        else:
            legend = Legend(items=legend_items)
            self.fig.add_layout(legend, 'center')

        self.fig.legend.click_policy = 'hide'    #- or 'mute'

    def create_zoomfig(self, viewer_cds, with_noise=True, with_coaddcam=True):
        #-----
        #- Zoom figure around mouse hover of main plot
        tooltips_zoomfig = [("wave","$x"),("flux","$y")]
        self.zoomfig = bk.figure(height=self.plot_height//2, width=self.plot_height//2,
            y_range=self.fig.y_range, x_range=(5000,5100),
            # output_backend="webgl",
            toolbar_location=None, tooltips=tooltips_zoomfig, tools=[])

        self.zoom_data_lines = list()
        self.zoom_noise_lines = list()
        for spec in viewer_cds.cds_spectra:
            self.zoom_data_lines.append(self.zoomfig.line('plotwave', 'plotflux', source=spec,
                line_color=self.colors[spec.name], line_width=1, line_alpha=self.alpha_discrete))
            if with_noise :
                self.zoom_noise_lines.append(self.zoomfig.line('plotwave', 'plotnoise', source=spec,
                                line_color=self.noise_colors[spec.name], line_width=1, line_alpha=self.alpha_discrete))
        if with_coaddcam :
            self.zoom_data_lines.append(self.zoomfig.line('plotwave', 'plotflux', source=viewer_cds.cds_coaddcam_spec, line_color=self.colors['coadd'], line_alpha=1))
            if with_noise :
                lx = self.zoomfig.line('plotwave', 'plotnoise', source=viewer_cds.cds_coaddcam_spec, line_color=self.noise_colors['coadd'], line_alpha=1)
                self.zoom_noise_lines.append(lx)

        if viewer_cds.cds_model is not None:
            lx = self.zoomfig.line('plotwave', 'plotflux', source=viewer_cds.cds_model, line_color='black')
        if viewer_cds.cds_othermodel is not None :
            lx = self.zoomfig.line('plotwave', 'plotflux', source=viewer_cds.cds_othermodel, line_color='black', line_dash='dashed')

        #- Callback to update zoom window x-range
        self.zoom_callback = CustomJS(
            args=dict(zoomfig=self.zoomfig,fig=self.fig),
            code="""
                zoomfig.x_range.start = cb_obj.x - 100;
                zoomfig.x_range.end = cb_obj.x + 100;
            """)

        self.fig.js_on_event(bokeh.events.MouseMove, self.zoom_callback)

    def create_imfig(self, spectra):
        #-----
        #- Targeting image
        npix_cutout = 256  # size of legacysurvey JPEG cutouts

        self.imfig = bk.figure(width=self.plot_height//2, height=self.plot_height//2,
                          x_range=(0, npix_cutout), y_range=(0, npix_cutout),
                          x_axis_location=None, y_axis_location=None,
                          output_backend="webgl",
                          toolbar_location=None, tools=[])
        self.imfig.min_border_left = 0
        self.imfig.min_border_right = 0
        self.imfig.min_border_top = 0
        self.imfig.min_border_bottom = 0

        self.imfig_urls = _viewer_urls(spectra)
        self.imfig_source = ColumnDataSource(data=dict(url=[self.imfig_urls[0][0]],
                                                       txt=[self.imfig_urls[0][2]]))
        imfig_img = self.imfig.image_url('url', source=self.imfig_source, x=1, y=1, w=256, h=256, anchor='bottom_left')
        imfig_txt = self.imfig.text(10, 256-30, text='txt', source=self.imfig_source,
                               text_color='yellow', text_font_size='8pt')

        #- This cross-hair is visible if the estimated PM correction is larger than 1 arcsec
        self.crosshair_source = ColumnDataSource(data=dict(xs=self.imfig_urls[0][3],
                                                           ys=self.imfig_urls[0][4]))
        pm_crosshair = self.imfig.multi_line('xs','ys', source=self.crosshair_source,
                                             line_width=1, line_color='white')

        #- Central cross-hair
        xs_center, ys_center = _cross_hair_points(npix_cutout//2, npix_cutout//2)
        central_crosshair = self.imfig.multi_line(xs_center, ys_center, line_width=1.5, line_color='yellow')


    def add_imfig_callback(self, viewer_widgets):
        #-----
        #- Targeting image callback
        # This has to be called once viewer_widgets are created => fct separated from create_imfig()
        self.imfig_callback = CustomJS(args = dict(
                                    urls = self.imfig_urls,
                                    ispectrumslider = viewer_widgets.ispectrumslider),
                                       code='''window.open(urls[ispectrumslider.value][1], "_blank");''')
        self.imfig.js_on_event('tap', self.imfig_callback)



    def add_spectral_lines(self, viewer_cds, figure='main', label_offset_top=100, label_offset_bottom=5):
        """Add spectral line markers to plot.

        Parameters
        ----------
        viewer_cds : array-like
            Viewer data.
        figure : {'main', 'zoom'}, optional
            Figure to add spectral lines to.
        label_offset_top : float, optional
            Offset in y-position for line labels with respect to top
            of the figure for emission lines.
        label_offset_bottom : float, optional
            Offset in y-position for line labels with respect to bottom
            of the figure for absorption lines.
        """

        if figure=='main' : bk_figure = self.fig
        elif figure=='zoom' : bk_figure = self.zoomfig
        else :
            raise ValueError("Unknown input figure type.")

        fig_height = bk_figure.plot_height
        if self.legend_outside_plot and figure=='main':
            label_offset_top += 10

        line_data = dict(viewer_cds.cds_spectral_lines.data)

        #- Labels y-position: default values
        y = list()
        for i in range(len(line_data['restwave'])):
            if line_data['emission'][i]:
                y.append(fig_height - label_offset_top)
            else:
                y.append(label_offset_bottom)

        #- Labels y-position: avoid overlaps
        restwave = np.asarray(line_data['restwave'])
        emission = np.asarray(line_data['emission'])
        for emission_flag in [True, False]:
            w, = np.where( (emission == emission_flag) )
            # Indices for the restwave array, sorted, and filtering emission_flag:
            sorted_indices = w[np.argsort(restwave[w])]
            for i in range(len(sorted_indices)-1):
                if (restwave[sorted_indices[i+1]] < restwave[sorted_indices[i]]+100):
                    # The following allows to have 4 levels of y-position offsets:
                    if emission_flag and (y[sorted_indices[i]]>=fig_height-label_offset_top-2*15):
                        y[sorted_indices[i+1]] = y[sorted_indices[i]] - 15
                    elif (not emission_flag) and (y[sorted_indices[i]]<=label_offset_bottom+2*15):
                        y[sorted_indices[i+1]] = y[sorted_indices[i]] + 15
        line_data['y'] = y

        #- Add vertical spans to figure
        if figure == 'main' :
            self.speclines = list()
            self.specline_labels = list()
        else :
            self.zoom_speclines = list()
            self.zoom_specline_labels = list()
        for w, y, name, emission, major in zip(
                line_data['plotwave'],
                line_data['y'],
                line_data['plotname'],
                line_data['emission'],
                line_data['major']
                ):
            color = 'blueviolet' if emission else 'green'
            visible = True if major else False

            s = Span(location=w, dimension='height', line_color=color,
                    line_alpha=1.0, line_dash='dashed', visible=visible)
            bk_figure.add_layout(s)

            lb = Label(x=w, y=y, x_units='data', y_units='screen',
                        text=name, text_color='gray', text_font_size="8pt",
                        x_offset=2, y_offset=0, visible=visible)
            bk_figure.add_layout(lb)

            if figure == 'main' :
                self.speclines.append(s)
                self.specline_labels.append(lb)
            else :
                self.zoom_speclines.append(s)
                self.zoom_specline_labels.append(lb)

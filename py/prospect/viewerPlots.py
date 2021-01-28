# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
====================
prospect.viewerPlots
====================

Class containing bokeh plots needed for the viewer

"""

import numpy as np
import bokeh.plotting as bk
from bokeh.models import ColumnDataSource, BoxAnnotation, Legend, Span, Label
import bokeh.layouts as bl
import bokeh.events

# TODO import _viewer_urls

class viewerPlots(object):
    """ 
    Encapsulates Bokeh plot-like objects that are part of prospect's GUI. 
    """
    # Copied from parts of viewer.plotspectra    

    def __init__(self):
        # Here I put, for convenience, all 'hardcoded' params...
        # All class members are not initialized here !!
        self.xmargin = 300.
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
        self.imfig = None

        # An attempt... (this is in the "targeting imaging" block in viewer.py)
        # (we could do similar for other "fig" members: Spacer by default)
        self.imfig = bl.Spacer(width=self.plot_height//2, height=self.plot_height//2)
        self.imfig_source = self.imfig_urls = None

    
    def create_mainfig(self, spectra, title, viewer_cds, sdss, with_noise=True, with_coaddcam=True):
        #-----
        #- Main figure
        #- Determine initial ymin, ymax, xmin, xmax
        self.ymin = self.ymax = self.xmax = 0
        self.xmin = 100000.
        if sdss:
            bands = ['coadd']
            self.ymin = np.nanmin(spectra.flux.value[0])
            self.ymax = np.nanmax(spectra.flux.value[0])
            self.xmin = np.min(spectra.spectral_axis.value)
            self.xmax = np.max(spectra.spectral_axis.value)
        else:
            bands = spectra.bands
            for i, band in enumerate(bands):
                if isinstance(spectra, SpectrumList):
                    sp_flux = spectra[i].flux.value[0]
                    sp_wave = spectra[i].spectral_axis.value
                else:
                    sp_flux = spectra.flux[band][0]
                    sp_wave = spectra.wave[band]
                self.ymin = min(self.ymin, np.nanmin(sp_flux))
                self.ymax = max(self.ymax, np.nanmax(sp_flux))
                self.xmin = min(self.xmin, np.min(sp_wave))
                self.xmax = max(self.xmax, np.max(sp_wave))
        self.xmin -= self.xmargin
        self.xmax += self.xmargin

        tools = 'pan,box_zoom,wheel_zoom,save'
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
                fill_alpha = self.alpha_overlapband if with_coaddcam else 0
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
        self.imfig = bk.figure(width=self.plot_height//2, height=self.plot_height//2,
                          x_range=(0, 256), y_range=(0, 256),
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
        imfig_txt = self.imfig.text(10, 256-30, text='txt', source=imfig_source,
                               text_color='yellow', text_font_size='8pt')
        # cross-hair
        self.imfig.multi_line([[129-15,129-5],[129+15,129+5],[129,129],[129,129]],
                         [[129,129],[129,129],[129-15,129-5],[129+5,129+15]], line_width=1, line_color='yellow')

    
    def add_spectral_lines(self, viewer_cds, figure='main', fig_height=None, label_offsets=[100, 5]):
        """
        label_offsets = [offset_absorption_lines, offset_emission_lines] : offsets in y-position
                        for line labels wrt top (resp. bottom) of the figure
        figure: 'main' or 'zoom' to flag if lines are added to self.fig or self.zoomfig
        """

        if figure=='main' : fig = self.fig
        elif figure=='zoom' : fig = self.zoomfig
        else :
            raise ValueError("Unknown input figure type.")

        if fig_height is None : fig_height = fig.plot_height

        line_data = dict(viewer_cds.cds_spectral_lines.data)
        y = list()
        for i in range(len(line_data['restwave'])):
            if i == 0:
                if line_data['emission'][i]:
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
        if figure == 'main' :
            self.speclines = list()
            self.specline_labels = list()
        else :
            self.zoom_speclines = list()
            self.zoom_specline_labels = list()
        for w, y, name, emission in zip(
                line_data['plotwave'],
                line_data['y'],
                line_data['plotname'],
                line_data['emission']
                ):
            if emission:
                color = 'blueviolet'
            else:
                color = 'green'

            s = Span(location=w, dimension='height', line_color=color,
                    line_alpha=1.0, line_dash='dashed', visible=False)

            lb = Label(x=w, y=y, x_units='data', y_units='screen',
                        text=name, text_color='gray', text_font_size="8pt",
                        x_offset=2, y_offset=0, visible=False)
                        
            if figure == 'main' :
                self.fig.add_layout(s)
                self.speclines.append(s)
                self.fig.add_layout(lb)
                self.specline_labels.append(lb)
            else :
                self.zoomfig.add_layout(s)
                self.zoom_speclines.append(s)
                self.zoomfig.add_layout(lb)
                self.zoom_specline_labels.append(lb)


        # TO THINK : name of class instation ? eg viewer_cds ? Au final c'est horrible viewer_cds.cds_spectra... 

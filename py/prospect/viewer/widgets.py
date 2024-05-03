# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
=======================
prospect.viewer.widgets
=======================

Class containing bokeh widgets needed for the viewer (except for VI widgets)

"""


import numpy as np

from bokeh.models import CustomJS, ColumnDataSource
from bokeh.models.widgets import (
    Slider, Button, TextInput, RadioButtonGroup, TableColumn,
    DataTable, CheckboxButtonGroup, CheckboxGroup, Select, Div)

from ..utilities import get_resources


def _metadata_table(table_keys, viewer_cds, table_width=500, shortcds_name='shortcds', selectable=False):
    """ Returns bokeh's (ColumnDataSource, DataTable) needed to display a set of metadata given by table_keys.

    """
    special_cell_width = { 'TARGETID':150, 'MORPHTYPE':70, 'SPECTYPE':70, 'SUBTYPE':60,
                         'Z':50, 'ZERR':50, 'Z_ERR':50, 'ZWARN':50, 'ZWARNING':50, 'DELTACHI2':70 }
    for x in viewer_cds.phot_bands:
        special_cell_width['mag_'+x] = 40
    special_cell_title = { 'DELTACHI2': 'Δχ2(N+1/N)' }

    table_columns = []
    cdsdata = dict()
    for key in table_keys:
        if key in special_cell_width.keys():
            cell_width = special_cell_width[key]
        else:
            cell_width = table_width//len(table_keys)
        if key in special_cell_title.keys():
            cell_title = special_cell_title[key]
        else:
            cell_title = key
        if ('mag_' in key) or ('EXPTIME' in key):
            cdsdata[key] = [ "{:.2f}".format(viewer_cds.cds_metadata.data[key][0]) ]
        elif 'CHI2' in key:
            cdsdata[key] = [ "{:.1f}".format(viewer_cds.cds_metadata.data[key][0]) ]
        elif key in ['Z', 'ZERR', 'Z_ERR']:
            cdsdata[key] = [ "{:.4f}".format(viewer_cds.cds_metadata.data[key][0]) ]
        else:
            cdsdata[key] = [ viewer_cds.cds_metadata.data[key][0] ]
        table_columns.append( TableColumn(field=key, title=cell_title, width=cell_width) )
    shortcds = ColumnDataSource(cdsdata, name=shortcds_name)
    # In order to be able to copy-paste the metadata in browser,
    #   the combination selectable=True, editable=True is needed:
    editable = True if selectable else False
    output_table = DataTable(source = shortcds, columns=table_columns,
                             index_position=None, selectable=selectable, editable=editable, width=table_width)
    output_table.height = 2 * output_table.row_height
    return (shortcds, output_table)


class ViewerWidgets(object):
    """
    Encapsulates Bokeh widgets, and related callbacks, that are part of prospect's GUI.
        Except for VI widgets
    """

    def __init__(self, plots, nspec):
        self.js_files = get_resources('js')
        self.navigation_button_width = 30
        self.z_button_width = 30
        self.plot_widget_width = (plots.plot_width+(plots.plot_height//2))//2 - 40 # used for widgets scaling

        #-----
        #- Ispectrumslider and smoothing widgets
        # Ispectrumslider's value controls which spectrum is displayed
        # These two widgets call update_plot(), later defined
        slider_end = nspec-1 if nspec > 1 else 0.5 # Slider cannot have start=end
        slidertitle = 'Spectrum number (0 to '+str(nspec-1)+')'
        self.ispectrumslider = Slider(start=0, end=slider_end, value=0, step=1, title=slidertitle)
        self.smootherslider = Slider(start=0, end=26, value=0, step=1.0, title='Gaussian Sigma Smooth')
        self.coaddcam_buttons = None
        self.model_select = None
        #- Small CDS to contain informations on widgets/plots status
        self.cds_widgetinfos = ColumnDataSource({'oii_save_xmin': [plots.fig.x_range.start],
                                                 'oii_save_xmax': [plots.fig.x_range.end],
                                                 'oii_save_nsmooth': [self.smootherslider.value],
                                                 'waveframe_active': [0],
                                                 'z_input_value': [0]
                                                })

    def add_navigation(self, nspec):
        #-----
        #- Navigation buttons
        self.prev_button = Button(label="<", width=self.navigation_button_width)
        self.next_button = Button(label=">", width=self.navigation_button_width)
        self.prev_callback = CustomJS(
            args=dict(ispectrumslider=self.ispectrumslider),
            code="""
            if(ispectrumslider.value>0 && ispectrumslider.end>=1) {
                ispectrumslider.value--
            }
            """)
        self.next_callback = CustomJS(
            args = dict(ispectrumslider=self.ispectrumslider, nspec=nspec),
            code = """
            if(ispectrumslider.value<nspec-1 && ispectrumslider.end>=1) {
                ispectrumslider.value++
            }
            """)
        self.prev_button.js_on_event('button_click', self.prev_callback)
        self.next_button.js_on_event('button_click', self.next_callback)
        #- Input spectrum number
        self.ispec_input = TextInput(value=str(self.ispectrumslider.value), width=50)
        self.ispec_input_callback = CustomJS(
            args = dict(ispec_input=self.ispec_input, ispectrumslider=self.ispectrumslider, nspec=nspec),
            code = """
            var i_spec = parseInt(ispec_input.value) ;
            if (Number.isInteger(i_spec) && i_spec>=0 && i_spec<nspec) {
                // Avoid recursive call
                if (i_spec != ispectrumslider.value) {
                    ispectrumslider.value = i_spec ;
                }
            }
            """)
        self.ispec_input.js_on_change('value', self.ispec_input_callback)

    def add_resetrange(self, viewer_cds, plots):
        #-----
        #- Axis reset button (superseeds the default bokeh "reset"
        self.reset_plotrange_button = Button(label="Reset X-Y range", button_type="default")
        reset_plotrange_code = self.js_files["adapt_plotrange.js"] + self.js_files["reset_plotrange.js"]
        self.reset_plotrange_callback = CustomJS(
            args = dict(fig = plots.fig,
                        xmin = plots.xmin,
                        xmax = plots.xmax,
                        spectra = viewer_cds.cds_spectra,
                        widgetinfos = self.cds_widgetinfos),
            code = reset_plotrange_code)
        self.reset_plotrange_button.js_on_event('button_click', self.reset_plotrange_callback)

    def add_redshift_widgets(self, z, viewer_cds, plots):
        ## TODO handle "z" (same issue as viewerplots TBD)

        #-----
        #- Redshift / wavelength scale widgets
        z1 = np.floor(z*100)/100
        dz = z-z1
        self.zslider = Slider(start=-0.1, end=5.0, value=z1, step=0.01, title='Redshift rough tuning')
        self.dzslider = Slider(start=0.0, end=0.0099, value=dz, step=0.0001, title='Redshift fine-tuning')
        self.zslider.format = "0[.]00" # default bokeh value, for record
        self.dzslider.format = "0[.]0000"
        self.z_input = TextInput(value="{:.4f}".format(z), title="Redshift value:")
        self.cds_widgetinfos.data['z_input_value'][0] = self.z_input.value

        #- Observer vs. Rest frame wavelengths
        self.waveframe_buttons = RadioButtonGroup(
            labels=["Obs", "Rest"], active=0)

        self.zslider_callback  = CustomJS(
            args=dict(zslider=self.zslider, dzslider=self.dzslider, z_input=self.z_input),
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

        self.dzslider_callback  = CustomJS(
            args=dict(zslider=self.zslider, dzslider=self.dzslider, z_input=self.z_input),
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

        self.zslider.js_on_change('value', self.zslider_callback)
        self.dzslider.js_on_change('value', self.dzslider_callback)

        self.z_minus_button = Button(label="<", width=self.z_button_width)
        self.z_plus_button = Button(label=">", width=self.z_button_width)
        self.z_minus_callback = CustomJS(
            args=dict(z_input=self.z_input),
            code="""
            var z = parseFloat(z_input.value)
            if(z >= -0.09) {
                z -= 0.01
                z_input.value = z.toFixed(4)
            }
            """)
        self.z_plus_callback = CustomJS(
            args=dict(z_input=self.z_input),
            code="""
            var z = parseFloat(z_input.value)
            if(z <= 4.99) {
                z += 0.01
                z_input.value = z.toFixed(4)
            }
            """)
        self.z_minus_button.js_on_event('button_click', self.z_minus_callback)
        self.z_plus_button.js_on_event('button_click', self.z_plus_callback)

        if 'Z' in viewer_cds.cds_metadata.data.keys():
            self.zreset_button = Button(label='Reset to z_pipe')
            self.zreset_callback = CustomJS(
                args=dict(z_input=self.z_input, metadata=viewer_cds.cds_metadata, ispectrumslider=self.ispectrumslider),
                code="""
                    var z = metadata.data['Z'][ispectrumslider.value]
                    z_input.value = z.toFixed(4)
                """)
            self.zreset_button.js_on_event('button_click', self.zreset_callback)
        else:
            self.zreset_button = Div(text="(z_pipe not available)")

        z_input_args = dict(spectra = viewer_cds.cds_spectra,
            coaddcam_spec = viewer_cds.cds_coaddcam_spec,
            model = viewer_cds.cds_model,
            othermodel = viewer_cds.cds_othermodel,
            metadata = viewer_cds.cds_metadata,
            widgetinfos = self.cds_widgetinfos,
            ispectrumslider = self.ispectrumslider,
            zslider = self.zslider,
            dzslider = self.dzslider,
            z_input = self.z_input,
            line_data = viewer_cds.cds_spectral_lines,
            lines = plots.speclines,
            line_labels = plots.specline_labels,
            zlines = plots.zoom_speclines,
            zline_labels = plots.zoom_specline_labels,
            overlap_waves = plots.overlap_waves,
            overlap_bands = plots.overlap_bands,
            fig = plots.fig)
        self.z_input_callback = CustomJS(
            args = z_input_args,
            code = self.js_files["shift_wave.js"] + self.js_files["change_redshift.js"]
        )
        self.z_input.js_on_change('value', self.z_input_callback)

        waveframe_args = z_input_args
        self.waveframe_callback = CustomJS(
            args = waveframe_args,
            code = self.js_files["shift_wave.js"] + self.js_files["change_waveframe.js"])
        self.waveframe_buttons.js_on_click(self.waveframe_callback)

    def add_oii_widgets(self, plots):
        #------
        #- Zoom on the OII doublet
        # TODO? optimize smoothing for autozoom (current value: 0)
        self.oii_zoom_button = Button(label="OII-zoom", button_type="default")
        self.oii_zoom_callback = CustomJS(
            args = dict(fig=plots.fig, smootherslider=self.smootherslider,
                        widgetinfos=self.cds_widgetinfos),
            code = """
            // Save previous setting (for the "Undo" button)
            widgetinfos.data['oii_save_xmin'][0] = fig.x_range.start
            widgetinfos.data['oii_save_xmax'][0] = fig.x_range.end
            widgetinfos.data['oii_save_nsmooth'][0] = smootherslider.value
            // Center on the middle of the redshifted OII doublet (vaccum)
            var central_wave = 3728.48;
            if (widgetinfos.data['waveframe_active'][0] == 0) {
                var z = parseFloat(widgetinfos.data['z_input_value'][0])
                central_wave *= (1+z)
            }
            fig.x_range.start = central_wave - 100
            fig.x_range.end = central_wave + 100
            // No smoothing (this implies a call to update_plot)
            smootherslider.value = 0
            """)
        self.oii_zoom_button.js_on_event('button_click', self.oii_zoom_callback)

        self.oii_undo_button = Button(label="Undo OII-zoom", button_type="default")
        self.oii_undo_callback = CustomJS(
            args = dict(fig=plots.fig, smootherslider=self.smootherslider, widgetinfos=self.cds_widgetinfos),
            code = """
            fig.x_range.start = widgetinfos.data['oii_save_xmin'][0]
            fig.x_range.end = widgetinfos.data['oii_save_xmax'][0]
            smootherslider.value = widgetinfos.data['oii_save_nsmooth'][0]
            """)
        self.oii_undo_button.js_on_event('button_click', self.oii_undo_callback)


    def add_coaddcam(self, plots):
        #-----
        #- Highlight individual-arm or camera-coadded spectra
        coaddcam_labels = ["Camera-coadded", "Single-arm"]
        self.coaddcam_buttons = RadioButtonGroup(labels=coaddcam_labels, active=0)
        self.coaddcam_callback = CustomJS(
            args = dict(coaddcam_buttons = self.coaddcam_buttons,
                        list_lines=[plots.data_lines, plots.noise_lines,
                                    plots.zoom_data_lines, plots.zoom_noise_lines],
                        alpha_discrete = plots.alpha_discrete,
                        overlap_bands = plots.overlap_bands,
                        alpha_overlapband = plots.alpha_overlapband),
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
        self.coaddcam_buttons.js_on_click(self.coaddcam_callback)


    def add_metadata_tables(self, viewer_cds, show_zcat=True,
                           top_metadata=['TARGETID', 'EXPID', 'COADD_NUMEXP', 'COADD_EXPTIME']):
        """ Display object-related informations
                top_metadata: metadata to be highlighted in table_a

            Note: "short" CDS, with a single row, are used to fill these bokeh tables.
            When changing object, js code modifies these short CDS so that tables are updated.
        """

        #- Sorted list of potential metadata:
        metadata_to_check = [ ('mag_'+x) for x in viewer_cds.phot_bands ]
        metadata_to_check += ['MORPHTYPE', 'TARGETID', 'HPXPIXEL', 'TILEID', 'COADD_NUMEXP', 'COADD_EXPTIME', 'COADD_NUMNIGHT',
                             'COADD_NUMTILE', 'NIGHT', 'EXPID', 'FIBER', 'CAMERA']
        table_keys = []
        for key in metadata_to_check:
            if key in viewer_cds.cds_metadata.data.keys():
                table_keys.append(key)
            if 'NUM_'+key in viewer_cds.cds_metadata.data.keys():
                for prefix in ['FIRST','LAST','NUM']:
                    table_keys.append(prefix+'_'+key)
                    if key in top_metadata:
                        top_metadata.append(prefix+'_'+key)

        #- Table a: "top metadata"
        table_a_keys = [ x for x in table_keys if x in top_metadata ]
        self.shortcds_table_a, self.table_a = _metadata_table(table_a_keys, viewer_cds, table_width=600,
                                                              shortcds_name='shortcds_table_a', selectable=True)
        #- Table b: Targeting information
        self.shortcds_table_b, self.table_b = _metadata_table(['Targeting masks'], viewer_cds, table_width=self.plot_widget_width,
                                                              shortcds_name='shortcds_table_b', selectable=True)
        #- Table(s) c/d : Other information (imaging, etc.)
        remaining_keys = [ x for x in table_keys if x not in top_metadata ]
        if len(remaining_keys) > 7:
            table_c_keys = remaining_keys[0:len(remaining_keys)//2]
            table_d_keys = remaining_keys[len(remaining_keys)//2:]
        else:
            table_c_keys = remaining_keys
            table_d_keys = None
        self.shortcds_table_c, self.table_c = _metadata_table(table_c_keys, viewer_cds, table_width=self.plot_widget_width,
                                                             shortcds_name='shortcds_table_c', selectable=False)
        if table_d_keys is None:
            self.shortcds_table_d, self.table_d = None, None
        else:
            self.shortcds_table_d, self.table_d = _metadata_table(table_d_keys, viewer_cds, table_width=self.plot_widget_width,
                                                                 shortcds_name='shortcds_table_d', selectable=False)

        #- Table z: redshift fitting information
        if show_zcat:
            if viewer_cds.dict_rrdetails is not None:  # "Detailled" info for Nth best fits
                fit_results = viewer_cds.dict_rrdetails
                # Case of DeltaChi2 : compute it from Chi2s
                #    The "DeltaChi2" in redrock files is between best fits for a given (spectype,subtype)
                #    Here we want to display DeltaChi2 independently of (spectype,subtype)
                #    Convention: DeltaChi2 = -1 for the last fit.
                chi2s = fit_results['CHI2'][0]
                full_deltachi2s = np.zeros(len(chi2s))-1
                full_deltachi2s[:-1] = chi2s[1:]-chi2s[:-1]
                cdsdata = dict(Nfit = np.arange(1,len(chi2s)+1),
                                SPECTYPE = fit_results['SPECTYPE'][0],  # [0:num_best_fits] (if we want to restrict... TODO?)
                                SUBTYPE = fit_results['SUBTYPE'][0],
                                Z = [ "{:.4f}".format(x) for x in fit_results['Z'][0] ],
                                ZERR = [ "{:.4f}".format(x) for x in fit_results['ZERR'][0] ],
                                ZWARN = fit_results['ZWARN'][0],
                                CHI2 = [ "{:.1f}".format(x) for x in fit_results['CHI2'][0] ],
                                DELTACHI2 = [ "{:.1f}".format(x) for x in full_deltachi2s ])
                self.shortcds_table_z = ColumnDataSource(cdsdata, name='shortcds_table_z')
                columns_table_z = [ TableColumn(field=x, title=t, width=w) for x,t,w in [ ('Nfit','Nfit',5), ('SPECTYPE','SPECTYPE',70), ('SUBTYPE','SUBTYPE',60), ('Z','Z',50) , ('ZERR','ZERR',50), ('ZWARN','ZWARN',50), ('DELTACHI2','Δχ2(N+1/N)',70)] ]
                self.table_z = DataTable(source=self.shortcds_table_z, columns=columns_table_z,
                                         selectable=False, index_position=None, width=self.plot_widget_width)
                self.table_z.height = 3 * self.table_z.row_height
            else :
                self.shortcds_table_z, self.table_z = _metadata_table(viewer_cds.zcat_keys, viewer_cds,
                                    table_width=self.plot_widget_width, shortcds_name='shortcds_table_z', selectable=False)
        else :
            self.table_z = Div(text="Not available ")
            self.shortcds_table_z = None


    def add_specline_toggles(self, viewer_cds, plots):
        #-----
        #- Toggle lines
        self.speclines_button_group = CheckboxButtonGroup(
                labels=["Emission lines", "Absorption lines"], active=[0, 1])
        self.majorline_checkbox = CheckboxGroup(
                labels=['Show only major lines'], active=[0])

        self.speclines_callback = CustomJS(
            args = dict(line_data = viewer_cds.cds_spectral_lines,
                        lines = plots.speclines,
                        line_labels = plots.specline_labels,
                        zlines = plots.zoom_speclines,
                        zline_labels = plots.zoom_specline_labels,
                        lines_button_group = self.speclines_button_group,
                        majorline_checkbox = self.majorline_checkbox),
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
        self.speclines_button_group.js_on_click(self.speclines_callback)
        self.majorline_checkbox.js_on_click(self.speclines_callback)


    def add_model_select(self, viewer_cds, num_approx_fits):
        #------
        #- Select secondary model to display
        model_options = []
        if viewer_cds.cds_model is not None:
            model_options = ['Best fit']
        if viewer_cds.cds_model_2ndfit is not None:
            model_options.append('2nd best fit')
        if num_approx_fits is not None:
            for i in range(1,1+num_approx_fits) :
                ith = 'th'
                if i==1 : ith='st'
                if i==2 : ith='nd'
                if i==3 : ith='rd'
                model_options.append(str(i)+ith+' fit (approx)')
        if viewer_cds.dict_std_templates is not None:
            std_template_labels = [x[5:] for x in viewer_cds.dict_std_templates.keys() if x[:5]=='wave_']
            for std_template in std_template_labels:
                model_options.append('STD '+std_template)
        self.model_select = Select(value=model_options[0],
                                   title="Other model (dashed curve):",
                                   options=model_options)
        model_select_code = self.js_files["interp_grid.js"] + self.js_files["smooth_data.js"] + self.js_files["select_model.js"]
        self.model_select_callback = CustomJS(
            args = dict(ispectrumslider = self.ispectrumslider,
                        model_select = self.model_select,
                        fit_templates = viewer_cds.dict_fit_templates,
                        cds_othermodel = viewer_cds.cds_othermodel,
                        cds_model_2ndfit = viewer_cds.cds_model_2ndfit,
                        cds_model = viewer_cds.cds_model,
                        rrdetails = viewer_cds.dict_rrdetails,
                        std_templates = viewer_cds.dict_std_templates,
                        median_spectra = viewer_cds.cds_median_spectra,
                        smootherslider = self.smootherslider,
                        z_input = self.z_input,
                        cds_metadata = viewer_cds.cds_metadata),
                        code = model_select_code)
        self.model_select.js_on_change('value', self.model_select_callback)


    def add_update_plot_callback(self, viewer_cds, plots, vi_widgets):
        #-----
        #- Main js code to update plots
        update_plot_code = (self.js_files["adapt_plotrange.js"] + self.js_files["interp_grid.js"] +
                            self.js_files["smooth_data.js"] + self.js_files["coadd_brz_cameras.js"] +
                            self.js_files["update_plot.js"])
        self.update_plot_callback = CustomJS(
            args = dict(
                spectra = viewer_cds.cds_spectra,
                coaddcam_spec = viewer_cds.cds_coaddcam_spec,
                model = viewer_cds.cds_model,
                othermodel = viewer_cds.cds_othermodel,
                model_2ndfit = viewer_cds.cds_model_2ndfit,
                metadata = viewer_cds.cds_metadata,
                rrdetails = viewer_cds.dict_rrdetails,
                shortcds_table_z = self.shortcds_table_z,
                shortcds_table_a = self.shortcds_table_a,
                shortcds_table_b = self.shortcds_table_b,
                shortcds_table_c = self.shortcds_table_c,
                shortcds_table_d = self.shortcds_table_d,
                ispectrumslider = self.ispectrumslider,
                ispec_input = self.ispec_input,
                smootherslider = self.smootherslider,
                z_input = self.z_input,
                widgetinfos = self.cds_widgetinfos,
                fig = plots.fig,
                xrange = [plots.xmin, plots.xmax],
                imfig_source = plots.imfig_source,
                crosshair_source = plots.crosshair_source,
                imfig_urls = plots.imfig_urls,
                model_select = self.model_select,
                vi_comment_input = vi_widgets.vi_comment_input,
                vi_std_comment_select = vi_widgets.vi_std_comment_select,
                vi_name_input = vi_widgets.vi_name_input,
                vi_quality_input = vi_widgets.vi_quality_input,
                vi_quality_labels = vi_widgets.vi_quality_labels,
                vi_issue_input = vi_widgets.vi_issue_input,
                vi_z_input = vi_widgets.vi_z_input,
                vi_category_select = vi_widgets.vi_category_select,
                vi_issue_slabels = vi_widgets.vi_issue_slabels
                ),
            code = update_plot_code
        )
        self.smootherslider.js_on_change('value', self.update_plot_callback)
        self.ispectrumslider.js_on_change('value', self.update_plot_callback)



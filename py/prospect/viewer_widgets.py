# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
=======================
prospect.viewer_widgets
=======================

Class containing bokeh widgets needed for the viewer (except for VI widgets)

"""


## TODO: when taking input class check isinstance() for all members needed in the function

import numpy as np

from bokeh.models import CustomJS, ColumnDataSource
from bokeh.models.widgets import (
    Slider, Button, TextInput, RadioButtonGroup, TableColumn,
    DataTable, CheckboxButtonGroup, CheckboxGroup, Select)

from .utilities import get_resources

class ViewerWidgets(object):
    """ 
    Encapsulates Bokeh widgets, and related callbacks, that are part of prospect's GUI.
        Except for VI widgets
    """
    
    def __init__(self, plots):
        self.js_files = get_resources('js')
        self.navigation_button_width = 30
        self.z_button_width = 30
        self.plot_widget_width = (plots.plot_width+(plots.plot_height//2))//2 - 40 # used for widgets scaling
    
        #-----
        #- Ifiberslider and smoothing widgets
        # Ifiberslider's value controls which spectrum is displayed
        # These two widgets call update_plot(), later defined
        slider_end = nspec-1 if nspec > 1 else 0.5 # Slider cannot have start=end
        self.ifiberslider = Slider(start=0, end=slider_end, value=0, step=1, title='Spectrum (of '+str(nspec)+')')
        self.smootherslider = Slider(start=0, end=26, value=0, step=1.0, title='Gaussian Sigma Smooth')
        self.coaddcam_buttons = None
        self.model_select = None


    def add_navigation(self):
        #-----
        #- Navigation buttons
        self.prev_button = Button(label="<", width=self.navigation_button_width)
        self.next_button = Button(label=">", width=self.navigation_button_width)
        self.prev_callback = CustomJS(
            args=dict(ifiberslider=self.ifiberslider),
            code="""
            if(ifiberslider.value>0 && ifiberslider.end>=1) {
                ifiberslider.value--
            }
            """)
        self.next_callback = CustomJS(
            args=dict(ifiberslider=self.ifiberslider, nspec=nspec),
            code="""
            if(ifiberslider.value<nspec-1 && ifiberslider.end>=1) {
                ifiberslider.value++
            }
            """)
        self.prev_button.js_on_event('button_click', self.prev_callback)
        self.next_button.js_on_event('button_click', self.next_callback)

    def add_resetrange(self, viewer_cds, plots):
        #-----
        #- Axis reset button (superseeds the default bokeh "reset"
        self.reset_plotrange_button = Button(label="Reset X-Y range", button_type="default")
        reset_plotrange_code = self.js_files["adapt_plotrange.js"] + self.js_files["reset_plotrange.js"]
        self.reset_plotrange_callback = CustomJS(args = dict(fig=plots.fig, xmin=plots.xmin, xmax=plots.xmax, spectra=viewer_cds.cds_spectra),
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
        self.dzslider.format = "0[.]0000"
        self.z_input = TextInput(value="{:.4f}".format(z), title="Redshift value:")

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

        self.zreset_button = Button(label='Reset to z_pipe')
        self.zreset_callback = CustomJS(
            args=dict(z_input=self.z_input, targetinfo=viewer_cds.cds_targetinfo, ifiberslider=self.ifiberslider),
            code="""
                var ifiber = ifiberslider.value
                var z = targetinfo.data['z'][ifiber]
                z_input.value = z.toFixed(4)
            """)
        self.zreset_button.js_on_event('button_click', self.zreset_callback)

        self.z_input_callback = CustomJS(
            args=dict(spectra = viewer_cds.cds_spectra,
                coaddcam_spec = viewer_cds.cds_coaddcam_spec,
                model = viewer_cds.cds_model,
                othermodel = viewer_cds.cds_othermodel,
                targetinfo = viewer_cds.cds_targetinfo,
                ifiberslider = self.ifiberslider,
                zslider = self.zslider,
                dzslider = self.dzslider,
                z_input = self.z_input,
                waveframe_buttons = self.waveframe_buttons,
                line_data = viewer_cds.cds_spectral_lines,
                lines = plots.speclines,
                line_labels = plots.specline_labels,
                zlines = plots.zoom_speclines,
                zline_labels = plots.zoom_specline_labels,
                overlap_waves = plots.overlap_waves,
                overlap_bands = plots.overlap_bands,
                fig = plots.fig
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
        self.z_input.js_on_change('value', self.z_input_callback)
        self.waveframe_buttons.js_on_click(self.z_input_callback)

        self.plotrange_callback = CustomJS(
            args = dict(
                z_input=self.z_input,
                waveframe_buttons=self.waveframe_buttons,
                fig=plots.fig,
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
        self.waveframe_buttons.js_on_click(self.plotrange_callback) # TODO: for record: is this related to waveframe bug? : 2 callbakcs for same click...


    def add_oii_widgets(self, plots):
        #------
        #- Zoom on the OII doublet TODO mv js code to other file
        # TODO: is there another trick than using a cds to pass the "oii_saveinfo" ?
        # TODO: optimize smoothing for autozoom (current value: 0)
        cds_oii_saveinfo = ColumnDataSource(
            {'xmin':[plots.fig.x_range.start], 'xmax':[plots.fig.x_range.end], 'nsmooth':[self.smootherslider.value]})
        self.oii_zoom_button = Button(label="OII-zoom", button_type="default")
        self.oii_zoom_callback = CustomJS(
            args = dict(z_input=self.z_input, fig=plots.fig, smootherslider=self.smootherslider,
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
        self.oii_zoom_button.js_on_event('button_click', self.oii_zoom_callback)

        self.oii_undo_button = Button(label="Undo", button_type="default")
        self.oii_undo_callback = CustomJS(
            args = dict(fig=plots.fig, smootherslider=self.smootherslider, cds_oii_saveinfo=cds_oii_saveinfo),
            code = """
            fig.x_range.start = cds_oii_saveinfo.data['xmin'][0]
            fig.x_range.end = cds_oii_saveinfo.data['xmax'][0]
            smootherslider.value = cds_oii_saveinfo.data['nsmooth'][0]
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
    
    
    def add_targetinfos(self, viewer_cds, sdss, show_zcat=True, template_dicts=None):
        #-----
        # Display object-related informations
        ## BYPASS DIV to be able to copy targetid...
        ## target_info_div = Div(text=cds_targetinfo.data['target_info'][0])
        tmp_dict = dict()
        tmp_dict['TARGETID'] = [ viewer_cds.cds_targetinfo.data['targetid'][0] ]
        tmp_dict['Target class'] = [ viewer_cds.cds_targetinfo.data['target_info'][0] ]
        targ_disp_cols = [ TableColumn(field='TARGETID', title='TARGETID', width=150),
                         TableColumn(field='Target class', title='Target class', width=250) ] # TODO tune width
        if sdss:
            phot_bands = ['u', 'g', 'r', 'i', 'z']
        else:
            phot_bands = ['G', 'R', 'Z', 'W1', 'W2']
        for band in phot_bands:
            tmp_dict['mag_'+band] = [ "{:.2f}".format(viewer_cds.cds_targetinfo.data['mag_'+band][0]) ]
            targ_disp_cols.append( TableColumn(field='mag_'+band, title='mag_'+band, width=40) )
        self.targ_disp_cds = ColumnDataSource(tmp_dict, name='targ_disp_cds')
        self.targ_display = DataTable(source = self.targ_disp_cds, columns=targ_disp_cols,index_position=None, selectable=True, editable=True) # width=...
        self.targ_display.height = 2 * self.targ_display.row_height
        if show_zcat is not None :
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
                tmp_dict = dict(SPECTYPE = [ viewer_cds.cds_targetinfo.data['spectype'][0] ],
                    SUBTYPE = [ viewer_cds.cds_targetinfo.data['subtype'][0] ],
                    Z = [ "{:.4f}".format(viewer_cds.cds_targetinfo.data['z'][0]) ],
                    ZERR = [ "{:.4f}".format(viewer_cds.cds_targetinfo.data['zerr'][0]) ],
                    ZWARN = [ viewer_cds.cds_targetinfo.data['zwarn'][0] ],
                    DeltaChi2 = [ "{:.1f}".format(viewer_cds.cds_targetinfo.data['deltachi2'][0]) ])
            self.zcat_disp_cds = ColumnDataSource(tmp_dict, name='zcat_disp_cds')
            zcat_disp_cols = [ TableColumn(field=x, title=t, width=w) for x,t,w in [ ('SPECTYPE','SPECTYPE',70), ('SUBTYPE','SUBTYPE',60), ('Z','Z',50) , ('ZERR','ZERR',50), ('ZWARN','ZWARN',50), ('DeltaChi2','Δχ2(N/N+1)',70)] ]
            if template_dicts is not None :
                zcat_disp_cols.insert(0, TableColumn(field='Nfit', title='Nfit', width=5))
            self.zcat_display = DataTable(source=self.zcat_disp_cds, columns=zcat_disp_cols, selectable=False, index_position=None, width=self.plot_widget_width)
            self.zcat_display.height = 2 * self.zcat_display.row_height
            if template_dicts is not None : self.zcat_display.height = 3 * self.zcat_display.row_height
        else :
            self.zcat_display = Div(text="Not available ")
            self.zcat_disp_cds = None


    def add_specline_toggles(self, viewer_cds, plots):
        #-----
        #- Toggle lines
        self.speclines_button_group = CheckboxButtonGroup(
                labels=["Emission lines", "Absorption lines"], active=[])
        self.majorline_checkbox = CheckboxGroup(
                labels=['Show only major lines'], active=[])

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
        self.lines_button_group.js_on_click(self.speclines_callback)
        self.majorline_checkbox.js_on_click(self.speclines_callback)


    def add_model_select(self, viewer_cds, template_dicts, num_approx_fits, with_full_2ndfit=True):
        #------
        #- Select secondary model to display
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
        self.model_select = Select(value=model_options[0], title="Other model (dashed curve):", options=model_options)
        model_select_code = self.js_files["interp_grid.js"] + self.js_files["smooth_data.js"] + self.js_files["select_model.js"]
        self.model_select_callback = CustomJS(
            args = dict(ifiberslider = self.ifiberslider,
                        model_select = self.model_select,
                        fit_templates=template_dicts[0],
                        cds_othermodel = viewer_cds.cds_othermodel,
                        cds_model_2ndfit = viewer_cds.cds_model_2ndfit,
                        cds_model = viewer_cds.cds_model,
                        fit_results=template_dicts[1],
                        std_templates=template_dicts[2],
                        median_spectra = viewer_cds.cds_median_spectra,
                        smootherslider = self.smootherslider,
                        z_input = self.z_input,
                        cds_targetinfo = viewer_cds.cds_targetinfo),
                        code = model_select_code)
        self.model_select.js_on_change('value', self.model_select_callback)


    def add_update_plot_callback(self, viewer_cds, plots, vi_widgets, template_dicts):
        #-----
        #- Main js code to update plots
        update_plot_code = (self.js_files["adapt_plotrange.js"] + self.js_files["interp_grid.js"] +
                            self.js_files["smooth_data.js"] + self.js_files["coadd_brz_cameras.js"] +
                            self.js_files["update_plot.js"])
        # TMP handling of template_dicts
        the_fit_results = None if template_dicts is None else template_dicts[1] # dirty
        self.update_plot_callback = CustomJS(
            args = dict(
                spectra = viewer_cds.cds_spectra,
                coaddcam_spec = viewer_cds.cds_coaddcam_spec,
                model = viewer_cds.cds_model,
                othermodel = viewer_cds.cds_othermodel,
                model_2ndfit = viewer_cds.cds_model_2ndfit,
                targetinfo = viewer_cds.cds_targetinfo,
                fit_results = the_fit_results,
                zcat_disp_cds = self.zcat_disp_cds,
                targ_disp_cds = self.targ_disp_cds,
                ifiberslider = self.ifiberslider,
                smootherslider = self.smootherslider,
                z_input = self.z_input,
                fig = plots.fig,
                imfig_source = plots.imfig_source,
                imfig_urls = plots.imfig_urls,
                model_select = self.model_select,
                vi_comment_input = vi_widgets.vi_comment_input,
                vi_std_comment_select = vi_widgets.vi_std_comment_select,
                vi_name_input = vi_widgets.vi_name_input,
                vi_class_input = vi_widgets.vi_class_input,
                vi_class_labels = vi_widgets.vi_class_labels,
                vi_issue_input = vi_widgets.vi_issue_input,
                vi_z_input = vi_widgets.vi_z_input,
                vi_category_select = vi_widgets.vi_category_select,
                vi_issue_slabels = vi_widgets.vi_issue_slabels
                ),
            code = update_plot_code
        )
        self.smootherslider.js_on_change('value', self.update_plot_callback)
        self.ifiberslider.js_on_change('value', self.update_plot_callback)



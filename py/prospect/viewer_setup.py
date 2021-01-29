# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
=====================
prospect.viewer_setup
=====================

Full bokeh setup for prospect

"""

import bokeh.layouts as bl
from bokeh.models import CustomJS, Tabs, Panel
from bokeh.models.widgets import Div
import bokeh.events

from .grid_thumbs import grid_thumbs

class ViewerSetup(object):
    '''
    Layout of plots and widgets for prospect's GUI
    '''
    # NB widget height / width are still partly hardcoded, but not arbitrary except for Spacers

    def __init__(self, plots, widgets, vi_widgets, with_vi_widgets=True):
        '''
        plots : :class:`ViewerPlots`
        widgets : :class:`ViewerWidgets`
        vi_widgets : :class:`ViewerVIWidgets`
        '''

        self.navigator = bl.row(
            bl.column(widgets.prev_button, width=widgets.navigation_button_width+15),
            bl.column(widgets.next_button, width=widgets.navigation_button_width+20),
            bl.column(widgets.ifiberslider, width=plots.plot_width+(plots.plot_height//2)-(60*len(vi_widgets.vi_class_labels)+2*widgets.navigation_button_width+35))
        )
        
        if with_vi_widgets :
            self.navigator.children.insert(1, bl.column(vi_widgets.vi_class_input, width=60*len(vi_widgets.vi_class_labels)) )
            if vi_widgets.vi_countdown_toggle is None :
                vi_header_block = bl.column( Div(text="VI optional indications :"), width=300 )
            else :
                vi_header_block = bl.row(
                    bl.column( Div(text="VI optional indications :"), width=300 ),
                    bl.column( vi_widgets.vi_countdown_toggle, width=200 )
                )
            self.vi_widget_set = bl.column(
                vi_header_block,
                bl.row(
                    bl.column(vi_widgets.vi_issue_input, width=150),
                    bl.column(vi_widgets.vi_z_input, width=150),
                    bl.column(vi_widgets.vi_category_select, width=150)
                ),
                bl.row(
                    bl.column(vi_widgets.vi_comment_input, width=300),
                    bl.column(vi_widgets.vi_std_comment_select, width=200),
                ),
                bl.row(
                    bl.column(vi_widgets.vi_name_input, width=200),
                    bl.column(vi_widgets.vi_filename_input, width=300)
                ),
                bl.column(vi_widgets.save_vi_button, width=100),
                bl.column(vi_widgets.vi_table),
                bl.row(
                    bl.column(vi_widgets.recover_vi_button, width=150),
                    bl.column(vi_widgets.clear_vi_button, width=150)
                ),
                background='#f5f5f0'
            )
            
        self.plot_widget_set = bl.column(
            bl.column( Div(text="Pipeline fit: ") ),
            bl.column(widgets.zcat_display, width=widgets.plot_widget_width),
            bl.row(
                bl.column(
                    bl.row(
                        bl.column(widgets.z_minus_button, width=widgets.z_button_width+15),
                        bl.column(widgets.zslider, width=widgets.plot_widget_width-2*widgets.z_button_width-135),
                        bl.column(widgets.z_plus_button, width=widgets.z_button_width)
                    ),
                    bl.row(
                        bl.column(widgets.dzslider, width=widgets.plot_widget_width-235),
                        bl.column(bl.Spacer(width=20)),
                        bl.column(widgets.zreset_button, width=100)
                    )
                ),
                bl.column(bl.Spacer(width=15)),
                bl.column(
                    bl.column(widgets.z_input, width=100),
                    bl.column(vi_widgets.z_tovi_button, width=100) # TODO if with_vi_widget
                ),
                background='#fff7e6'
            ),
            bl.column(widgets.smootherslider, width=widgets.plot_widget_width),
            bl.row(
                bl.column(widgets.speclines_button_group, width=200),
                bl.column(bl.Spacer(width=30)),
                bl.column(widgets.majorline_checkbox, width=120)
            )
        )
        if widgets.coaddcam_buttons is not None :
            waveframe_block = bl.row(
                                bl.column(widgets.coaddcam_buttons, width=200),
                                bl.column(bl.Spacer(width=30)),
                                bl.column(widgets.waveframe_buttons, width=120)
                              )
        else :
            waveframe_block = bl.column(widgets.waveframe_buttons, width=120)
        self.plot_widget_set.children.append(waveframe_block)
        if widgets.model_select is not None :
            self.plot_widget_set.children.insert(3, bl.column(widgets.model_select, width=200))
        
        if with_vi_widgets :
            self.plot_widget_set.children.append( bl.column(bl.Spacer(height=30)) )
            self.plot_widget_set.children.append( bl.column(vi_widgets.vi_guideline_div, width=widgets.plot_widget_width) )
            self.full_widget_set = bl.row(
                self.vi_widget_set,
                bl.column(bl.Spacer(width=40)),
                self.plot_widget_set
            )
        else : self.full_widget_set = self.plot_widget_set

        self.main_bokehsetup = bl.column(
            bl.row(plots.fig, bl.column(plots.imfig, plots.zoomfig), bl.Spacer(width=20)),
            bl.row(
                bl.column(widgets.targ_display, width=600), # plot_width - 200
                bl.column(bl.Spacer(width=20)),
                bl.column(widgets.reset_plotrange_button, width = 120),
                bl.column(bl.Spacer(width=80)),
                bl.column(widgets.oii_zoom_button, width=80),
                bl.column(widgets.oii_undo_button, width=50),
            ),
            self.navigator,
            self.full_widget_set,
            sizing_mode='stretch_width'
        )

        self.full_viewer = self.main_bokehsetup


    def add_thumb_tab(self, spectra, plots, widgets, nspec):
        self.ncols_grid = 5 # TODO un-hardcode
        self.miniplot_width = ( plots.plot_width + (plots.plot_height//2) ) // self.ncols_grid

        self.full_viewer = Tabs()
        titles = None # TODO define
        self.thumb_grid = grid_thumbs(spectra, self.miniplot_width,
                x_range=(plots.xmin,plots.xmax),
                ncols_grid=self.ncols_grid, titles=titles)
        tab1 = Panel(child = self.main_bokehsetup, title='Main viewer')
        tab2 = Panel(child = self.thumb_grid, title='Gallery')
        self.full_viewer.tabs=[ tab1, tab2 ]

        # Dirty trick : callback functions on thumbs need to be defined AFTER the full_viewer is implemented
        # Otherwise, at least one issue = no toolbar anymore for main fig. (apparently due to ifiberslider in callback args)
        for i_spec in range(nspec) :
            self.thumb_callback = CustomJS(args=dict(full_viewer=self.full_viewer, i_spec=i_spec, ifiberslider=widgets.ifiberslider), code="""
            full_viewer.active = 0
             ifiberslider.value = i_spec
            """)
            (self.thumb_grid.children[i_spec][0]).js_on_event(bokeh.events.DoubleTap, self.thumb_callback)



class ThumbGrid(object):
    ## Standalone grid of simple thumbs (to make lightweight pages showing all spectra)
    
    def __init__(self, spectra, plots, title):
        self.ncols_grid = 5 # TODO un-hardcode
        titles = None # TODO define
        self.miniplot_width = (plots.plot_width + (plots.plot_height//2) ) // ncols_grid
        self.thumb_grid = grid_thumbs(spectra, miniplot_width, x_range=(plots.xmin,plots.xmax), ncols_grid=ncols_grid, titles=titles)
        self.thumb_viewer = bl.column(
            bl.column( Div(text=
                           " <h3> Thumbnail gallery for DESI spectra in "+title+" </h3>" +
                           " <p> Click <a href='specviewer_"+title+".html'>here</a> to access the spectral viewer corresponding to these spectra. </p>"
                          ), width=plots.plot_width ),
            bl.column( thumb_grid )
        )


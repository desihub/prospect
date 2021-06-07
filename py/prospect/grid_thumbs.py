# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
====================
prospect.grid_thumbs
====================

Grid of thumbnail plots displaying a set of spectra

"""

import numpy as np
import astropy.convolution

import bokeh.plotting as bk
import bokeh.layouts as bl

from .coaddcam import coaddcam_prospect

# TODOS:
# - structure: class??
# - improve thumbs with infos, plot quality

def grid_thumbs(spectra, thumb_width, x_range=(3400,10000), thumb_height=None, resamp_factor=15, ncols_grid=5, titles=None) :
    '''
    Create a bokeh gridplot of thumbnail pictures from spectra
    - coadd arms
    - smooth+resample to reduce size of embedded CDS, according to resamp_factor
    - titles : optional list of titles for each thumb

    TODO: Not tested on Spectrum1D objects.
    '''

    if thumb_height is None : thumb_height = thumb_width//2
    if titles is not None : assert len(titles) == spectra.num_spectra()
    thumb_wave, thumb_flux, dummy = coaddcam_prospect(spectra)
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


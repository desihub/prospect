# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
===================
prospect.coaddcam
===================

Python routines which translate the simple camera-coaddition algorithms
used in prospect webpages (js code in js/interp_grid.js and js/coadd_brz_cameras.js).
"""

import numpy as np
from math import floor

def index_dichotomy(point, grid):
    """Find nearest index in `grid`, left from `point`; use dichotomy method.

    Translated from js/interp_grid.js.

    Parameters
    ----------
    point : :class:`float`
        Value to find in `grid`.
    grid : array-like
        Values to search.

    Returns
    -------
    :class:`int`
        Nearest index.
    """
    if point < grid[0]:
        return 0
    if point > grid[-1]:
        return len(grid) - 2
    i_left = 0
    i_center = 0
    i_right = len(grid) - 1
    while i_right - i_left != 1:
        i_center = i_left + floor((i_right - i_left)/2)
        if point >= grid[i_center]:
            i_left = i_center
        else:
            i_right = i_center
    return i_left


def interp_grid(xval, xarr, yarr):
    """Basic linear interpolation of [`xarr`, `yarr`] on point `xval`.

    Translated from js/interp_grid.js.

    Parameters
    ----------
    xval : :class:`xval`
        Interpolate y-value at this point.
    xarr, yarr : array-like
        X, Y data.

    Returns
    -------
    :class:`float`
        The y-value corresponding to `xval`.
    """
    index = index_dichotomy(xval, xarr)
    a = (yarr[index+1] - yarr[index])/(xarr[index+1] - xarr[index])
    b = yarr[index]-a*xarr[index]
    yval = a*xval+b
    return yval


def coadd_brz_cameras(wave_in, flux_in, noise_in) :
    """Camera-coadd *brz* spectra.

    Translated from js/coadd_brz_cameras.js.

    Parameters
    ----------
    wave_in : array-like
        Set of three wavelength arrays corresponding to *brz*.
    flux_in : array-like
        Set of three flux arrays corresponding to *brz*.
    noise_in : array-like
        Noise arrays for weighting.

    Returns
    -------
    :func:`tuple`
        The coadded wavelength solution, flux and noise.

    Notes
    -----
    * In case of a missing arm (or 2), the data are just concatenated.
    * Need to handle case of no noise.
    """
    wave_out = []
    flux_out = []
    noise_out = []

    # Special case of missing arm, handled separately:
    if len(wave_in)<=2:
        for i_cam in range(len(wave_in)):
            for i in range(len(wave_in[i_cam])):
                wave_out.append(wave_in[i_cam][i])
                flux_out.append(flux_in[i_cam][i])
                noise_out.append(noise_in[i_cam][i])
        return (np.asarray(wave_out), np.asarray(flux_out), np.asarray(noise_out))

    # Find b,r,z ordering in input arrays
    wave_start = [wave_in[0][0], wave_in[1][0], wave_in[2][0]]
    i_b = wave_start.index(np.amin(wave_start))
    i_z = wave_start.index(np.amax(wave_start))
    i_r = 1
    for i in [0,1,2] :
        if ( (i_b != i) and (i_z != i) ) : i_r = i

    margin = 20
    for i in range(len(wave_in[i_b])) : # b
        if (wave_in[i_b][i] < wave_in[i_b][-1] - margin) :
            wave_out.append(wave_in[i_b][i])
            flux_out.append(flux_in[i_b][i])
            noise_out.append(noise_in[i_b][i])
    the_lim = wave_out[-1]
    for i in range(len(wave_in[i_r])) : # r
        if ( (wave_in[i_r][i] < wave_in[i_r][-1] - margin) and (wave_in[i_r][i] > the_lim)) :
            wave_out.append(wave_in[i_r][i])
            flux_out.append(flux_in[i_r][i])
            noise_out.append(noise_in[i_r][i])
    the_lim = wave_out[-1]
    for i in range(len(wave_in[i_z])) : # z
        if (wave_in[i_z][i] > the_lim) :
            wave_out.append(wave_in[i_z][i])
            flux_out.append(flux_in[i_z][i])
            noise_out.append(noise_in[i_z][i])
    for i in range(len(wave_out)) : # combine in overlapping regions
        b1 = -1
        b2 = -1
        if ( (wave_out[i] > wave_in[i_r][0]) and (wave_out[i] < wave_in[i_b][-1]) ) : # br
            b1 = 0
            b2 = 1
        if ( (wave_out[i] > wave_in[i_z][0]) and (wave_out[i] < wave_in[i_r][-1]) ) : # rz
            b1 = 1
            b2 = 2
        if (b1 != -1) :
            phi1 = interp_grid(wave_out[i], wave_in[b1], flux_in[b1])
            noise1 = interp_grid(wave_out[i], wave_in[b1], noise_in[b1])
            phi2 = interp_grid(wave_out[i], wave_in[b2], flux_in[b2])
            noise2 = interp_grid(wave_out[i], wave_in[b2], noise_in[b2])
            if ( noise1 > 0 and noise2 > 0 ) :
                iv1 = 1/(noise1*noise1)
                iv2 = 1/(noise2*noise2)
                iv = iv1+iv2
                noise_out[i] = 1/np.sqrt(iv)
                flux_out[i] = (iv1*phi1+iv2*phi2)/iv
    return (np.asarray(wave_out), np.asarray(flux_out), np.asarray(noise_out))


def coaddcam_prospect(spectra):
    """ Camera-coaddition of *brz* bands in a set of DESI spectra.

    This is essentially a wrapper to :func:`~prospect.coaddcam.coadd_brz_cameras`.

    Parameters
    ----------
    spectra: :class:`~desispec.spectra.Spectra` or similar

    Returns
    -------
    :func:`tuple` (wave, flux, ivar), where
        wave : 1D[nwave] array of wavelengths
        flux : 2D[nspec, nwave] array of flux densities
        ivar : 2D[nspec, nwave] array of inverse variances of `flux`

    Raises
    ------
    RuntimeError
        If spectra.bands does not contain 'b', 'r', 'z', and spectra.bands is not 'brz'
    """

    if np.any([ band in spectra.bands for band in ['b','r','z'] ]) :

        for i_spec in range(spectra.num_spectra()) :
            wave_in = [ spectra.wave[b] for b in spectra.bands ]
            flux_in = [ spectra.flux[b][i_spec] for b in spectra.bands ]
            noise_in = []
            for b in spectra.bands :
                the_noise = np.zeros(len(spectra.ivar[b][i_spec]))
                w, = np.where( (spectra.ivar[b][i_spec] > 0) )
                the_noise[w] = 1/np.sqrt(spectra.ivar[b][i_spec][w])
                noise_in.append(the_noise)
            the_wave, the_flux, the_noise = coadd_brz_cameras(wave_in, flux_in, noise_in)
            the_ivar = np.zeros(len(the_noise))
            w, = np.where( (the_noise > 0) )
            the_ivar[w] = 1/the_noise[w]**2
            if i_spec == 0 :
                wave_out = np.asarray(the_wave)
                flux_out = np.zeros((spectra.num_spectra(),len(wave_out)),
                                    dtype = list(spectra.flux.values())[0].dtype)
                ivar_out = np.zeros((spectra.num_spectra(),len(wave_out)),
                                    dtype = list(spectra.ivar.values())[0].dtype)
            flux_out[i_spec,:] = the_flux
            ivar_out[i_spec,:] = the_ivar

    elif spectra.bands == ['brz'] :
        wave_out = spectra.wave['brz']
        flux_out = spectra.flux['brz']
        ivar_out = spectra.ivar['brz']

    else :
        raise RuntimeError("Set of bands for spectra not supported.")

    return (wave_out, flux_out, ivar_out)

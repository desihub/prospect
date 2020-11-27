
# EA - Oct 2019 (Temporary / preliminary)
# desispec.coaddition.coadd_cameras() unsatisfying at least since
# 1) don't want to coadd over exposures / 2) cannot assume waves are aligned over arms (r/z mismatch seen in datachallenge)

# In addition we need a python implementation of the same algorithm as used in js/coadd_brz_cameras.js

import numpy as np
import math

from desispec.interpolation import resample_flux

def index_dichotomy(point, grid) :
    """
    Translated from js/interp_grid.js
    Find nearest index in grid, left from point; use dichotomy method
    """
    if ( point < grid[0] ) : return 0
    if ( point > grid[-1] ) : return len(grid)-2
    i_left = 0
    i_center = 0
    i_right = len(grid)-1
    while ( i_right - i_left != 1) :
        i_center = i_left + math.floor((i_right-i_left)/2)
        if ( point >= grid[i_center] ) :
            i_left = i_center
        else :
            i_right = i_center
    return i_left


def interp_grid(xval, xarr, yarr) :
    """
    Translated from js/interp_grid.js
    Basic linear interpolation of [xarr,yarr] on point xval
    """
    index = index_dichotomy(xval, xarr)
    a = (yarr[index+1] - yarr[index])/(xarr[index+1] - xarr[index])
    b = yarr[index]-a*xarr[index]
    yval = a*xval+b
    return yval


def coadd_brz_cameras(wave_in, flux_in, noise_in) :
    """
    Translated from js/coadd_brz_cameras.js
    Camera-coadd brz spectra.
        each "_in" must have 3 entries (brz)
        TODO handle case of no noise
    """

    # Find b,r,z ordering in input arrays
    wave_start = [wave_in[0][0], wave_in[1][0], wave_in[2][0]]
    i_b = wave_start.index(np.amin(wave_start))
    i_z = wave_start.index(np.amax(wave_start))
    i_r = 1
    for i in [0,1,2] :
        if ( (i_b != i) and (i_z != i) ) : i_r = i

    wave_out = []
    flux_out = []
    noise_out = []
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


def coaddcam_prospect(spectra) :
    """
    Wrapper to coadd_brz_cameras.
    Same input/output as mycoaddcam()
    """
    
    if np.all([ band in spectra.bands for band in ['b','r','z'] ]) :
    
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
                flux_out = np.zeros((spectra.num_spectra(),len(wave_out)),dtype=spectra.flux['b'].dtype)
                ivar_out = np.zeros((spectra.num_spectra(),len(wave_out)),dtype=spectra.ivar['b'].dtype)
            flux_out[i_spec,:] = the_flux
            ivar_out[i_spec,:] = the_ivar

    elif spectra.bands == ['brz'] :
        wave_out = spectra.wave['brz']
        flux_out = spectra.flux['brz']
        ivar_out = spectra.ivar['brz']
    else :
        raise RuntimeError("Set of bands for spectra not supported.")
 
    return (wave_out, flux_out, ivar_out)


def mycoaddcam(spectra) :
    """"
    Merges brz spectra into a single (wave,flux)
      takes into account noise and mis-matched wavelengths over the 3 arms
    Currently assumes b r z bands and two overlap regions
    """
    
    if np.all([ band in spectra.bands for band in ['b','r','z'] ]) :
    
        # Define (arbitrarily) wavelength grid
        margin = 20 # Angstrom. Avoids using edge-of-band at overlap regions
        wave = spectra.wave['b'].copy()
        wave = wave[ (wave<np.max(wave)-margin) ]
        tolerance = 0.0001
        length_bands = {'b' : wave.size}
        w_bands = {'b' : np.arange(wave.size)}
        for band in ['r','z'] :
            if band=='z' : w_bands[band], = np.where( spectra.wave[band]>wave[-1]+tolerance )
            else : w_bands[band], = np.where( (spectra.wave[band]>wave[-1]+tolerance)
                                      & (spectra.wave[band]<np.max(spectra.wave[band])-margin) )
            wave=np.append(wave,spectra.wave[band][w_bands[band]])
            length_bands[band] = w_bands[band].size

        nwave = wave.size
        nspec = spectra.num_spectra()
        flux = np.zeros((nspec,nwave),dtype=spectra.flux['b'].dtype)
        ivar = np.zeros((nspec,nwave),dtype=spectra.ivar['b'].dtype)

        # Flux in non-overlapping waves
        i = 0
        for band in ['b', 'r', 'z'] :
            flux[:,i:i+length_bands[band]] = spectra.flux[band][:,w_bands[band]]
            ivar[:,i:i+length_bands[band]] = spectra.ivar[band][:,w_bands[band]]
            i += length_bands[band]

        # Overlapping regions
        overlaps = ['br','rz']
        for the_overlap in overlaps :
            b1, b2 = the_overlap[0], the_overlap[1]
            w_overlap, = np.where( (wave > spectra.wave[b2][0]) & (wave < spectra.wave[b1][-1]) )
            assert (w_overlap.size > 0)
            lambd_over = wave[w_overlap]
            for ispec in range(nspec) :
                phi1, ivar1 = resample_flux(lambd_over, spectra.wave[b1], spectra.flux[b1][ispec,:], ivar=spectra.ivar[b1][ispec,:])
                phi2, ivar2 = resample_flux(lambd_over, spectra.wave[b2], spectra.flux[b2][ispec,:], ivar=spectra.ivar[b2][ispec,:])
                ivar[ispec,w_overlap] = ivar1+ivar2
                w_ok = np.where( ivar[ispec,w_overlap] > 0)
                flux[ispec,w_overlap] = (phi1+phi2)/2
                flux[ispec,w_overlap][w_ok] = (ivar1[w_ok]*phi1[w_ok] + ivar2[w_ok]*phi2[w_ok])/ivar[ispec,w_overlap][w_ok]
    
    elif spectra.bands == ['brz'] :
        wave = spectra.wave['brz']
        flux = spectra.flux['brz']
        ivar = spectra.ivar['brz']
    else :
        raise RuntimeError("mycoaddcam: set of bands for spectra not supported")
    
    
    return (wave, flux, ivar)

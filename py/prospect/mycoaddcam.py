
# EA - Oct 2019 (Temporary / preliminary)
# desispec.coaddition.coadd_cameras() unsatisfying at least since
# 1) don't want to coadd over exposures / 2) cannot assume waves are aligned over arms (r/z mismatch seen in datachallenge)

import numpy as np

from desispec.interpolation import resample_flux

def mycoaddcam(spectra) :
    """"
    Merges brz spectra into a single (wave,flux)
      takes into account noise and mis-matched wavelengths over the 3 arms
    Currently assumes b r z bands and two overlap regions
    """
    
    assert np.all([ band in spectra.wave.keys() for band in ['b','r','z'] ]) 
    
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
    
    return (wave, flux, ivar)

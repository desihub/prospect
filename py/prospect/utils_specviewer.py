# -*- coding: utf-8 -*-

"""
Utility functions for prospect
"""

import numpy as np
import astropy.io.fits
from astropy.table import Table, vstack
import scipy.ndimage.filters


import matplotlib
matplotlib.use('Agg') # No DISPLAY
import matplotlib.pyplot as plt

from desiutil.log import get_logger
from prospect import mycoaddcam

_vi_flags = [
    # Definition of VI flags
    # Replaces former list viflags = ["Yes","No","Maybe","LowSNR","Bad"]
    # shortlabels for "issue" flags must be a unique single-letter identifier
    {"label" : "4", "type" : "class", "description" : "Confident classification, two or more secure features"},
    {"label" : "3", "type" : "class", "description" : "Probable classification, at least one secure feature + continuum; or many weak features"},
    {"label" : "2", "type" : "class", "description" : "Possible classification, one strong emission feature, but not sure what it is"},
    {"label" : "1", "type" : "class", "description" : "Unlikely classification, one or some unidentified features"},
    {"label" : "0", "type" : "class", "description" : "Nothing there"},
    {"label" : "Bad redshift fit", "shortlabel" : "R", "type" : "issue", "description" : "Misidentification of redshift by pipeline fitter"},
    {"label" : "Bad spectype fit", "shortlabel" : "C", "type" : "issue", "description" : "Misidentification of spectral category (eg. star vs QSO...)"},
    {"label" : "Bad spectrum", "shortlabel" : "S", "type" : "issue", "description" : "Bad spectrum (eg. cosmic/skyline residuals...)"}
]

_vi_file_fields = [
    # Contents of VI files: [ "field name (in VI file header)", "associated variable in cds_targetinfo", default dtype in VI file]
    # Ordered list
    ["TargetID", "targetid", "i4"],
    ["ExpID", "expid", "i4"],
    ["Spec version", "spec_version", "i4"], # TODO define
    ["Redrock version", "redrock_version", "i4"], # TODO define
    ["Redrock spectype", "spectype", "S6"],
    ["Redrock z", "z", "f4"],
    ["VI scanner", "VI_scanner", "S6"],
    ["VI class", "VI_class_flag", "i4"],
    ["VI issue", "VI_issue_flag", "S6"],
    ["VI z", "VI_z", "f4"],
    ["VI spectype", "VI_spectype", "S6"],
    ["VI comment", "VI_comment", "S6"]
]

_vi_spectypes =[
    # List of spectral types to fill in VI categories
    # in principle, it should match somehow redrock spectypes...
    "STAR",
    "GALAXY",
    "QSO"
]

def read_vi(vifile) :
    '''
    Read visual inspection file (ASCII/CSV or FITS according to file extension)
    Return full VI catalog, in Table format
    '''
    vi_records = [x[0] for x in _vi_file_fields]

    if (vifile[-5:] != ".fits" and vifile[-4:] not in [".fit",".fts",".csv"]) :
        raise RuntimeError("wrong file extension")
    if vifile[-4:] == ".csv" :
        vi_info = Table.read(vifile,format='ascii.csv', names=vi_records)
    else :
        vi_info = astropy.io.fits.getdata(vifile,1)
        if [(x in vi_info.names) for x in vi_records]!=[1 for x in vi_records] :
            raise RuntimeError("wrong records in vi fits file")
        vi_info = Table(vi_info)

    return vi_info


def match_vi_targets(vifile, targetlist) :
    '''
    Returns list of VIs matching the list of targetids
    For a given target, several VI entries can be available
    '''
    vi_info = read_vi(vifile)
    vicatalog=[ [] for i in range(len(targetlist)) ]
    for itarget,targetnum in enumerate(targetlist) :
        w,=np.where( (vi_info['targetid'] == targetnum) )
        if len(w)>0 : vicatalog[itarget] = vi_info[w]
    return vicatalog


def convert_vi_tofits(vifile_in, overwrite=True) :
    log = get_logger()
    if vifile_in[-4:] != ".csv" : raise RuntimeError("wrong file extension")
    vi_info = read_vi(vifile_in)
    vifile_out=vifile_in.replace(".csv",".fits")
    vi_info.write(vifile_out, format='fits', overwrite=overwrite)
    log.info("Created fits file : "+vifile_out+" ("+str(len(vi_info))+" entries).")
    

def initialize_master_vi(mastervifile, overwrite=False) :
    '''
    Create "master" VI file with no entry
    '''
    log = get_logger()
    vi_records = [x[0] for x in _vi_file_fields]
    vi_dtype = [x[2] for x in _vi_file_fields]
    vi_info = Table(names=vi_records, dtype=tuple(vi_dtype))
    vi_info.write(mastervifile, format='fits', overwrite=overwrite)
    log.info("Initialized VI file : "+mastervifile+" (0 entry)")
    

def merge_vi(mastervifile, newvifile) :
    '''
    Merge a new VI file to the "master" VI file
    The master file is overwritten.
    '''
    log = get_logger()
    mastervi = read_vi(mastervifile)
    newvi = read_vi(newvifile)
    mergedvi = vstack([mastervi,newvi], join_type='exact')
    mergedvi.write(mastervifile, format='fits', overwrite=True)
    log.info("Updated master VI file : "+mastervifile+" (now "+str(len(mergedvi))+" entries).")


def match_zcat_to_spectra(zcat_in, spectra) :
    '''
    zcat_in : astropy Table from redshift fitter
    creates a new astropy Table whose rows match the targetids of input spectra
    also returns the corresponding list of indices
    '''
    zcat_out = Table(dtype=zcat_in.dtype)
    index_list = list()
    for i_spec in range(spectra.num_spectra()) :
        ww, = np.where((zcat_in['TARGETID'] == spectra.fibermap['TARGETID'][i_spec]))
        if len(ww)<1 : raise RuntimeError("zcat table cannot match spectra.")
        zcat_out.add_row(zcat_in[ww[0]])
        index_list.append(ww[0])
    return (zcat_out, index_list)


def get_y_minmax(pmin, pmax, data, ispec) :
    '''
    Utility, from plotframe
    '''
    dx = np.sort(data[np.isfinite(data)])
    if len(dx)==0 : return (0,0)
    imin = int(np.floor(pmin*len(dx)))
    imax = int(np.floor(pmax*len(dx)))
    if (imax >= len(dx)) : imax = len(dx)-1
    return (dx[imin],dx[imax])


def miniplot_spectrum(spectra, i_spec, model=None, saveplot=None, smoothing=-1, coaddcam=True) :
    '''
    Matplotlib version of plotspectra, to plot a given spectrum
    Pieces of code were copy-pasted from plotspectra()
    Smoothing option : simple gaussian filtering
    '''
    data=[]
    if coaddcam is True :
        wave, flux, ivar = mycoaddcam.mycoaddcam(spectra)
        bad = ( ivar == 0.0 )
        flux[bad] = np.nan
        thedat = dict(
                band = 'coadd',
                wave = wave,
                flux = flux[i_spec]
                )
        data.append(thedat)
    else :
        for band in spectra.bands :
            #- Set masked bins to NaN so that Bokeh won't plot them
            bad = (spectra.ivar[band] == 0.0) | (spectra.mask[band] != 0)
            spectra.flux[band][bad] = np.nan
            thedat=dict(
                band = band,
                wave = spectra.wave[band].copy(),
                flux = spectra.flux[band][i_spec].copy()
                )
            data.append(thedat)
    
    if model is not None:
        mwave, mflux = model
        mwave = mwave.copy() # in case of
        mflux = mflux[i_spec].copy()

    # Gaussian smoothing
    if smoothing > 0 :
        ymin,ymax=0,0
        for spec in data :
            spec['wave'] = spec['wave']
            spec['flux'] = scipy.ndimage.filters.gaussian_filter1d(spec['flux'], sigma=smoothing, mode='nearest')
            tmpmin,tmpmax=get_y_minmax(0.01, 0.99, spec['flux'],i_spec)
            ymin=np.nanmin((tmpmin,ymin))
            ymax=np.nanmax((tmpmax,ymax))
        if model is not None :
            mwave = mwave[int(smoothing):-int(smoothing)]
            mflux = scipy.ndimage.filters.gaussian_filter1d(mflux, sigma=smoothing, mode='nearest')[int(smoothing):-int(smoothing)]
    
    colors = dict(b='#1f77b4', r='#d62728', z='maroon', coadd='#d62728')
    # for visibility, do not plot near-edge of band data (noise is high there):
    waverange = dict(b=[3500,5800], r=[5800,7600], z=[7600,9900], coadd=[3500,9900])
    for spec in data :
        band = spec['band']
        w, = np.where( (spec['wave']>=waverange[band][0]) & (spec['wave']<=waverange[band][1]) )
        plt.plot(spec['wave'][w],spec['flux'][w],c=colors[band])
    if model is not None :
        plt.plot(mwave, mflux, c='k')
    # No label to save space
    if smoothing > 0 :
        ymin = ymin*1.4 if (ymin<0) else ymin*0.6
        ymax = ymax*1.4
        plt.ylim((ymin,ymax))
    # TODO : include some infos on plot
    
    if saveplot is not None : plt.savefig(saveplot, dpi=50) # default dpi=100, TODO tune dpi
    else : print("No plot saved?") # TODO saveplot kwd optional or not ?
    plt.clf()
        
    return
    
    

# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
===================
prospect.viewer.cds
===================

Class containing all bokeh's ColumnDataSource objects needed in viewer.py

"""

import numpy as np
from pkg_resources import resource_filename

import bokeh.plotting as bk
from bokeh.models import ColumnDataSource

_specutils_imported = True
try:
    from specutils import Spectrum1D, SpectrumList
except ImportError:
    _specutils_imported = False

_desitarget_imported = True
try:
    from desitarget.targetmask import desi_mask
    from desitarget.cmx.cmx_targetmask import cmx_mask
    from desitarget.sv1.sv1_targetmask import desi_mask as sv1_desi_mask
    from desitarget.sv1.sv1_targetmask import bgs_mask as sv1_bgs_mask
except ImportError:
    _desitarget_imported = False

from ..mycoaddcam import coaddcam_prospect


def _airtovac(w):
    """Convert air wavelengths to vacuum wavelengths. Don't convert less than 2000 Å.

    Parameters
    ----------
    w : :class:`float`
        Wavelength [Å] of the line in air.

    Returns
    -------
    :class:`float`
        Wavelength [Å] of the line in vacuum.
    """
    if w < 2000.0:
        return w;
    vac = w
    for iter in range(2):
        sigma2 = (1.0e4/vac)*(1.0e4/vac)
        fact = 1.0 + 5.792105e-2/(238.0185 - sigma2) + 1.67917e-3/(57.362 - sigma2)
        vac = w*fact
    return vac


class ViewerCDS(object):
    """
    Encapsulates Bokeh ColumnDataSource objects to be passed to js callback functions.
    """
    
    def __init__(self):
        self.cds_spectra = None
        self.cds_median_spectra = None
        self.cds_coaddcam_spec = None
        self.cds_model = None
        self.cds_model_2ndfit = None
        self.cds_othermodel = None
        self.cds_metadata = None
    
    def load_spectra(self, spectra, with_noise=True):
        """ Creates column data source for observed spectra """
        
        self.cds_spectra = list()
        is_desispec = False
        if _specutils_imported and isinstance(spectra, SpectrumList):
            s = spectra
            bands = spectra.bands
        elif _specutils_imported and isinstance(spectra, Spectrum1D):
            s = [spectra]
            bands = ['coadd']
        else : # Assume desispec Spectra obj
            is_desispec = True
            s = spectra
            bands = spectra.bands
        
        for j, band in enumerate(bands):
            input_wave = s.wave[band] if is_desispec else s[j].spectral_axis.value
            input_nspec = spectra.num_spectra() if is_desispec else s[j].flux.shape[0]
            cdsdata = dict(
                origwave = input_wave.copy(),
                plotwave = input_wave.copy(),
                )
            for i in range(input_nspec):
                key = 'origflux'+str(i)
                input_flux = spectra.flux[band][i] if is_desispec else s[j].flux.value[i, :]
                cdsdata[key] = input_flux.copy()
                if with_noise :
                    key = 'orignoise'+str(i)
                    input_ivar = spectra.ivar[band][i] if is_desispec else s[j].uncertainty.array[i, :]
                    noise = np.zeros(len(input_ivar))
                    w, = np.where( (input_ivar > 0) )
                    noise[w] = 1/np.sqrt(input_ivar[w])
                    cdsdata[key] = noise
            cdsdata['plotflux'] = cdsdata['origflux0']
            if with_noise : 
                cdsdata['plotnoise'] = cdsdata['orignoise0']
            self.cds_spectra.append( ColumnDataSource(cdsdata, name=band) )
    
    def compute_median_spectra(self, spectra):
        """ Stores the median value for each spectrum into CDS.
            Simple concatenation of all values from different bands.
        """
                
        cdsdata = dict(median=[])
        for i in range(spectra.num_spectra()):
            flux_array = np.concatenate( tuple([spectra.flux[band][i] for band in spectra.bands]) )
            w, = np.where( ~np.isnan(flux_array) )
            if len(w)==0 :
                cdsdata['median'].append(1)
            else :
                cdsdata['median'].append(np.median(flux_array[w]))

        self.cds_median_spectra = ColumnDataSource(cdsdata)
        
    def init_coaddcam_spec(self, spectra, with_noise=True):
        """ Creates column data source for camera-coadded observed spectra
            Do NOT store all coadded spectra in CDS obj, to reduce size of html files
            Except for the first spectrum, coaddition is done later in javascript
        """
        
        coadd_wave, coadd_flux, coadd_ivar = coaddcam_prospect(spectra)
        cds_coaddcam_data = dict(
            origwave = coadd_wave.copy(),
            plotwave = coadd_wave.copy(),
            plotflux = coadd_flux[0,:].copy(),
            plotnoise = np.ones(len(coadd_wave))
        )
        if with_noise :
            w, = np.where( (coadd_ivar[0,:] > 0) )
            cds_coaddcam_data['plotnoise'][w] = 1/np.sqrt(coadd_ivar[0,:][w])
        self.cds_coaddcam_spec = ColumnDataSource(cds_coaddcam_data)

    def init_model(self, model, second_fit=False):
        """ Creates a CDS for model spectrum """
        
        mwave, mflux = model
        cdsdata = dict(
            origwave = mwave.copy(),
            plotwave = mwave.copy(),
            plotflux = np.zeros(len(mwave)),
        )
        for i in range(len(mflux)):
            key = 'origflux'+str(i)
            cdsdata[key] = mflux[i]
        cdsdata['plotflux'] = cdsdata['origflux0']
        
        if second_fit:
            self.cds_model_2ndfit = ColumnDataSource(cdsdata)
        else:
            self.cds_model = ColumnDataSource(cdsdata)

    def init_othermodel(self, zcatalog):
        """ Initialize CDS for the 'other model' curve, from the best fit """
        self.cds_othermodel = ColumnDataSource({
            'plotwave' : self.cds_model.data['plotwave'],
            'origwave' : self.cds_model.data['origwave'],
            'origflux' : self.cds_model.data['origflux0'],
            'plotflux' : self.cds_model.data['origflux0'],
            'zref' : zcatalog['Z'][0]+np.zeros(len(self.cds_model.data['origflux0'])) # Track z reference in model
        })
    
    
    def load_metadata(self, spectra, mask_type=None, zcatalog=None, username=" ", survey='DESI'):
        """ Creates column data source for target-related metadata,
            from fibermap, zcatalog and VI files 
        """
    
        if survey == 'DESI':
            nspec = spectra.num_spectra()
            self.zcat_keys = ['Z', 'SPECTYPE', 'SUBTYPE', 'ZERR', 'ZWARN', 'DELTACHI2']
            self.phot_bands = ['G','R','Z', 'W1', 'W2']
            assert _desitarget_imported
            supported_masks = {
                'SV1_DESI_TARGET': sv1_desi_mask,
                'SV1_BGS_TARGET': sv1_bgs_mask,
                'DESI_TARGET': desi_mask,
                'CMX_TARGET': cmx_mask
                }
        elif survey == 'SDSS':
            nspec = spectra.flux.shape[0]
            self.zcat_keys = ['Z', 'CLASS', 'SUBCLASS', 'Z_ERR', 'ZWARNING', 'RCHI2DIFF']
            self.phot_bands = ['u', 'g', 'r', 'i', 'z']
            supported_masks = ['PRIMTARGET', 'SECTARGET',
                                'BOSS_TARGET1', 'BOSS_TARGET2',
                                'ANCILLARY_TARGET1', 'ANCILLARY_TARGET2',
                                'EBOSS_TARGET0', 'EBOSS_TARGET1', 'EBOSS_TARGET2']
        else:
            raise ValueError('Wrong survey')
        
        self.cds_metadata = ColumnDataSource()
        # TODO: everywhere: handle case of no metadata: 
        #     Do something or not ?? it's better if we can avoid, need to check
        
        #- Generic metadata
        if survey == 'DESI':
            fibermap_keys = ['EXPID', 'NIGHT', 'TILEID'] #, 'MORPHTYPE']
            for fm_key in fibermap_keys:
                if fm_key in spectra.fibermap.keys():
                    self.cds_metadata.add(spectra.fibermap[fm_key], name=fm_key) 
               # else:
                    # ? self.cds_targetinfo.add(['-1' for i in range(nspec)], name=cds_key)
            #- Special case for targetids: No int64 in js !!
            self.cds_metadata.add([str(x) for x in spectra.fibermap['TARGETID']], name='TARGETID')
        elif survey == 'SDSS':
            self.cds_metadata.add([str(x.tolist()) for x in spectra.meta['plugmap']['OBJID']], name='TARGETID')
        
        #- Photometry
        for i, bandname in enumerate(self.phot_bands) :
            # TODO: define what todo if FLX_BAND not in fibermap
            if survey == 'SDSS':
                mag = spectra.meta['plugmap']['MAG'][:, i]
            else :
                mag = np.zeros(nspec)
                flux = spectra.fibermap['FLUX_'+bandname]
                extinction = np.ones(len(flux))
                if ('MW_TRANSMISSION_'+bandname) in spectra.fibermap.keys() :
                    extinction = spectra.fibermap['MW_TRANSMISSION_'+bandname]
                w, = np.where( (flux>0) & (extinction>0) )
                mag[w] = -2.5*np.log10(flux[w]/extinction[w])+22.5
            self.cds_metadata.add(mag, name='mag_'+bandname)

        #- Targeting masks
        if mask_type is not None:
            if survey == 'DESI':
                if mask_type not in spectra.fibermap.keys():
                    raise ValueError("mask_type is not in spectra.fibermap: "+mask_type)
                mask_used = supported_masks[mask_type]
                target_bits = spectra.fibermap[mask_type]
                target_info = [ ' '.join(mask_used.names(x)) for x in target_bits ]
            elif survey == 'SDSS':
                assert mask_type in supported_masks
                target_info = [ mask_type + ' (DUMMY)' for x in spectra.meta['plugmap'] ] # placeholder
            self.cds_metadata.add(target_info, name='target_info')

        #- Software versions
        #- TODO : get redrock version (from zcatalog...)
        if survey == 'SDSS':
            desispec_version = 'SDSS'
        else :
            desispec_specversion = "0"
            for xx,yy in spectra.meta.items() :
                if yy=="desispec" : desispec_specversion = spectra.meta[xx.replace('NAM','VER')]
        self.cds_metadata.add([desispec_specversion for i in range(nspec)], name='spec_version')
        self.cds_metadata.add(np.zeros(nspec)-1, name='redrock_version')
        self.cds_metadata.add(np.zeros(nspec)-1, name='template_version')

        #- Redshift fit
        if zcatalog is not None:
            for zcat_key in self.zcat_keys:
                if 'TYPE' in zcat_key or 'CLASS' in zcat_key:
                    data = zcatalog[zcat_key].astype('U{0:d}'.format(zcatalog[zcat_key].dtype.itemsize))
                else :
                    data = zcatalog[zcat_key]
                self.cds_metadata.add(data, name=zcat_key)
        
        #- VI informations
        vi_values = {
            'VI_scanner': username,
            'VI_class_flag': "-1",
            'VI_issue_flag': "",
            'VI_z': "",
            'VI_spectype': "",
            'VI_comment': ""
            } # this could go to utilities
        for vi_key, vi_val in vi_values.items():
            self.cds_metadata.add([vi_val for i in range(nspec)], name=vi_key)
            
    
    def load_spectral_lines(self, z=0):    
    
        line_data = dict(
            restwave = [],
            plotwave = [],
            name = [],
            longname = [],
            plotname = [],
            emission = [],
            major = [],
            #y = []
        )
        for line_category in ('emission', 'absorption'):
            # encoding=utf-8 is needed to read greek letters
            line_array = np.genfromtxt(resource_filename('prospect', "data/{0}_lines.txt".format(line_category)),
                                       delimiter=",",
                                       dtype=[("name", "|U20"),
                                              ("longname", "|U20"),
                                              ("wavelength", float),
                                              ("vacuum", bool),
                                              ("major", bool)],
                                        encoding='utf-8')
            vacuum_wavelengths = line_array['wavelength']
            w, = np.where(line_array['vacuum']==False)
            vacuum_wavelengths[w] = np.array([_airtovac(wave) for wave in line_array['wavelength'][w]])
            line_data['restwave'].extend(vacuum_wavelengths)
            line_data['plotwave'].extend(vacuum_wavelengths * (1+z))
            line_data['name'].extend(line_array['name'])
            line_data['longname'].extend(line_array['longname'])
            line_data['plotname'].extend(line_array['name'])
            emission_flag = True if line_category=='emission' else False
            line_data['emission'].extend([emission_flag for row in line_array])
            line_data['major'].extend(line_array['major'])

        self.cds_spectral_lines = ColumnDataSource(line_data)





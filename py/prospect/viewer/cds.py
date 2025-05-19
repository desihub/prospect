# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
===================
prospect.viewer.cds
===================

Class containing all bokeh's ColumnDataSource objects needed in viewer.py

"""
import importlib.resources
import numpy as np
from astropy.io import fits
from astropy.table import Table

import bokeh.plotting as bk
from bokeh.models import ColumnDataSource

_specutils_imported = True
try:
    from specutils import Spectrum1D, SpectrumList
except ImportError:
    _specutils_imported = False

_desispec_imported = True
try:
    from desispec.interpolation import resample_flux
except ImportError:
    _desispec_imported = False

from ..coaddcam import coaddcam_prospect
from ..utilities import supported_desitarget_masks, vi_file_fields, load_redrock_templates


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
        return w
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
        self.cds_spectral_lines = None
        self.dict_fit_templates = None  # Special case: not a CDS
        self.dict_std_templates = None  # Special case: not a CDS
        self.dict_rrdetails = None  # Special case: not a CDS

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
        is_desispec = False
        if _specutils_imported and isinstance(spectra, SpectrumList):
            s = spectra
            bands = spectra.bands
            nspec = spectra[0].flux.shape[0]
        elif _specutils_imported and isinstance(spectra, Spectrum1D):
            s = [spectra]
            bands = ['coadd']
            nspec = spectra.flux.shape[0]
        else : # Assume desispec Spectra obj
            is_desispec = True
            s = spectra
            bands = spectra.bands
            nspec = spectra.num_spectra()
        for i in range(nspec):
            if is_desispec:
                flux_array = np.concatenate( tuple([s.flux[band][i] for band in bands]) )
            else:
                flux_array = np.concatenate( tuple([s[j].flux[i, :].value for j, band in enumerate(bands)]) )
            w, = np.where( ~np.isnan(flux_array) )
            if len(w)==0 :
                cdsdata['median'].append(1.0)
            else :
                cdsdata['median'].append(np.median(flux_array[w]).tolist())

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


    def load_fit_templates(self, template_dir=None, nbpts_templates=4000, zcat_header=None):
        """ Create dict for spectral templates used in Redrock fits.
            These are used to recompute Redrock's Nth best-fit spectra on-the-fly
            in javascript.
            Templates are resampled in order to limit the size of html pages (and the
            browser's CPU usage).
            This resampling is dictated by parameter nbpts_templates.
            zcat_header is header from Redrock output with TEMNAMnn/TEMVERnn keywords
            indicating the version of the templates used at the time of the fit.
        """
        assert _desispec_imported # for resample_flux
        rr_templts = load_redrock_templates(template_dir=template_dir, zcat_header=zcat_header)
        self.dict_fit_templates = dict()
        for key,templt in rr_templts.items():
            fulltype_key = "_".join(key)   # merge redrock's (TYPE, SUBTYPE)
            wave_array = np.linspace(templt.wave[0], templt.wave[-1], num=nbpts_templates)
            flux_array = np.zeros(( templt.flux.shape[0],len(wave_array) ))
            for i in range(templt.flux.shape[0]):
                flux_array[i,:] = resample_flux(wave_array, templt.wave, templt.flux[i,:])
            self.dict_fit_templates["wave_"+fulltype_key] = wave_array
            self.dict_fit_templates["flux_"+fulltype_key] = flux_array


    def load_std_templates(self, std_template_file=None):
        """ Load a dict of "standard" templates.
            The std template file is `data/std_templates.fits`.
            It was created from `../scripts/prospect_std_templates.py`.
        """
        self.dict_std_templates = dict()
        if std_template_file is None:
            std_template_file = importlib.resources.files('prospect').joinpath("data", "std_templates.fits")
        hdul = fits.open(std_template_file)
        nhdu = len(hdul)
        hdul.close()
        for i in range(1, nhdu):
            t = Table.read(std_template_file, hdu=i)
            for key in t.keys():
                #- check table column name:
                if key[:5] not in ['wave_', 'flux_']:
                    raise ValueError('STD template file: wrong column name ('+key+')')
                #- check wavelength array is regularly, linearly binned (with absolute tolerance 0.01 AA):
                if key[:5]=='wave_':
                    waves = np.array(t[key])
                    delta_waves = waves[1:] - waves[:-1]
                    if not np.allclose(delta_waves, delta_waves[0], atol=0.01, rtol=1.e-10):
                        raise ValueError('STD template file: found irregular wavelength binning ('+key+')')
                self.dict_std_templates[key] = np.array(t[key])
        #- initialize cds_othermodel, if this was not done yet:
        if self.cds_othermodel is None:
            key_zero = list(self.dict_std_templates.keys())[0][5:]
            self.cds_othermodel = ColumnDataSource({
                'plotwave' : self.dict_std_templates['wave_'+key_zero],
                'origwave' : self.dict_std_templates['wave_'+key_zero],
                'origflux' : self.dict_std_templates['flux_'+key_zero],
                'plotflux' : self.dict_std_templates['flux_'+key_zero],
                'zref' : np.zeros(len(self.dict_std_templates['flux_'+key_zero]))  # std templates have z=0
            })


    def load_rrdetails(self, redrock_cat):
        """ Create dict for detailled redrock outputs.
            Used to recompute redrock's Nth best fit spectra on-the-fly in javascript,
            and display them in a table.
        """
        self.dict_rrdetails = dict()
        for key in redrock_cat.keys() :
            self.dict_rrdetails[key] = np.asarray(redrock_cat[key])
        self.dict_rrdetails['Nfit'] = redrock_cat['Z'].shape[1]


    def load_metadata(self, spectra, mask_type=None, zcatalog=None, survey='DESI'):
        """ Creates column data source for target-related metadata,
            from fibermap, zcatalog and VI files
        """

        if survey == 'DESI':
            nspec = spectra.num_spectra()
            # Optional metadata:
            fibermap_keys = ['HPXPIXEL', 'MORPHTYPE', 'CAMERA',
                             'COADD_NUMEXP', 'COADD_EXPTIME',
                             'COADD_NUMNIGHT', 'COADD_NUMTILE']
            # Optional metadata, will check matching FIRST/LAST/NUM keys in fibermap:
            special_fm_keys = ['FIBER', 'NIGHT', 'EXPID', 'TILEID']
            # Mandatory keys if zcatalog is set:
            self.zcat_keys = ['Z', 'SPECTYPE', 'SUBTYPE', 'ZERR', 'ZWARN', 'DELTACHI2']
            # Mandatory metadata:
            self.phot_bands = ['G','R','Z', 'W1', 'W2']
            supported_masks = supported_desitarget_masks
            # Galactic extinction coefficients:
            # - Wise bands from https://github.com/dstndstn/tractor/blob/master/tractor/sfd.py
            # - Other bands from desiutil.dust (updated coefficients Apr 2021,
            #   matching https://desi.lbl.gov/trac/wiki/ImagingStandardBandpass)
            R_extinction = {'W1':0.184, 'W2':0.113, 'W3':0.0241, 'W4':0.00910,
                            'G_N':3.258, 'R_N':2.176, 'Z_N':1.199,
                            'G_S':3.212, 'R_S':2.164, 'Z_S':1.211}
        elif survey == 'SDSS':
            nspec = spectra.flux.shape[0]
            # Mandatory keys if zcatalog is set:
            self.zcat_keys = ['Z', 'CLASS', 'SUBCLASS', 'Z_ERR', 'ZWARNING', 'RCHI2DIFF']
            # Mandatory metadata:
            self.phot_bands = ['u', 'g', 'r', 'i', 'z']
            supported_masks = ['PRIMTARGET', 'SECTARGET',
                                'BOSS_TARGET1', 'BOSS_TARGET2',
                                'ANCILLARY_TARGET1', 'ANCILLARY_TARGET2',
                                'EBOSS_TARGET0', 'EBOSS_TARGET1', 'EBOSS_TARGET2']
        else:
            raise ValueError('Wrong survey')

        self.cds_metadata = ColumnDataSource()

        #- Generic metadata
        if survey == 'DESI':
            #- Special case for targetids: No int64 in js !!
            self.cds_metadata.add([str(x) for x in spectra.fibermap['TARGETID']], name='TARGETID')
            #- "Special" keys: check for FIRST/LAST/NUM
            for fm_key in special_fm_keys:
                use_first_last_num = False
                if all([ (x+fm_key in spectra.fibermap.keys()) for x in ['FIRST_','LAST_','NUM_'] ]):
                    if np.any(spectra.fibermap['NUM_'+fm_key] > 1) : # if NUM==1, use fm_key only
                        use_first_last_num = True
                        self.cds_metadata.add(spectra.fibermap['FIRST_'+fm_key], name='FIRST_'+fm_key)
                        self.cds_metadata.add(spectra.fibermap['LAST_'+fm_key], name='LAST_'+fm_key)
                        self.cds_metadata.add(spectra.fibermap['NUM_'+fm_key], name='NUM_'+fm_key)
                if (not use_first_last_num) and fm_key in spectra.fibermap.keys():
                    # Do not load placeholder metadata:
                    if not (np.all(spectra.fibermap[fm_key]==0) or np.all(spectra.fibermap[fm_key]==-1)):
                        self.cds_metadata.add(spectra.fibermap[fm_key], name=fm_key)
            #- "Normal" keys
            for fm_key in fibermap_keys:
                # Arbitrary choice:
                if fm_key == 'COADD_NUMEXP' and 'NUM_EXPID' in self.cds_metadata.data.keys():
                    continue
                if fm_key == 'COADD_NUMNIGHT' and 'NUM_NIGHT' in self.cds_metadata.data.keys():
                    continue
                if fm_key == 'COADD_NUMTILE' and 'NUM_TILEID' in self.cds_metadata.data.keys():
                    continue
                if fm_key in spectra.fibermap.keys():
                    if not (np.all(spectra.fibermap[fm_key]==0) or np.all(spectra.fibermap[fm_key]==-1)):
                        self.cds_metadata.add(spectra.fibermap[fm_key], name=fm_key)
        elif survey == 'SDSS':
            #- Set 'TARGETID' name to OBJID for convenience
            self.cds_metadata.add([str(x.tolist()) for x in spectra.meta['plugmap']['OBJID']], name='TARGETID')

        #- Photometry
        for i, bandname in enumerate(self.phot_bands) :
            if survey == 'SDSS':
                mag = spectra.meta['plugmap']['MAG'][:, i]
            else :
                mag = np.zeros(nspec)
                flux = spectra.fibermap['FLUX_'+bandname]
                extinction = np.ones(len(flux))
                if ('MW_TRANSMISSION_'+bandname) in spectra.fibermap.keys():
                    extinction = spectra.fibermap['MW_TRANSMISSION_'+bandname]
                elif ('EBV' in spectra.fibermap.keys()) and (bandname.upper() in ['W1','W2','W3','W4']):
                    extinction = 10**(- R_extinction[bandname.upper()] * spectra.fibermap['EBV'])
                elif all(x in spectra.fibermap.keys() for x in ['EBV','PHOTSYS']) and (bandname.upper() in ['G','R','Z']):
                    for photsys in ['N', 'S']:
                        wphot, = np.where(spectra.fibermap['PHOTSYS'] == photsys)
                        a_band = R_extinction[bandname.upper()+"_"+photsys] * spectra.fibermap['EBV'][wphot]
                        extinction[wphot] = 10**(-a_band / 2.5)
                w, = np.where( (flux>0) & (extinction>0) )
                mag[w] = -2.5*np.log10(flux[w]/extinction[w])+22.5
            self.cds_metadata.add(mag, name='mag_'+bandname)

        #- Targeting masks
        if mask_type is not None:
            if survey == 'DESI':
                if mask_type not in spectra.fibermap.keys():
                    mask_candidates = [x for x in spectra.fibermap.keys() if '_TARGET' in x]
                    raise ValueError(mask_type+" is not in spectra.fibermap.\n Hints of available masks: "+(' '.join(mask_candidates)))
                mask_used = supported_masks[mask_type]
                target_bits = spectra.fibermap[mask_type]
                target_info = [ ' '.join(mask_used.names(x)) for x in target_bits ]
            elif survey == 'SDSS':
                assert mask_type in supported_masks
                target_info = [ mask_type + ' (DUMMY)' for x in spectra.meta['plugmap'] ] # placeholder
            self.cds_metadata.add(target_info, name='Targeting masks')

        #- Software versions
        #- TODO : get template version (from zcatalog...)
        if survey == 'SDSS':
            spec_version = 'SDSS'
        else :
            spec_version = '0'
            for xx,yy in spectra.meta.items() :
                if yy=="desispec" : spec_version = spectra.meta[xx.replace('NAM','VER')]
        self.cds_metadata.add([spec_version for i in range(nspec)], name='spec_version')
        redrock_version = ["-1" for i in range(nspec)]
        if zcatalog is not None:
            if 'RRVER' in zcatalog.keys(): redrock_version = zcatalog['RRVER'].data
        self.cds_metadata.add(redrock_version, name='redrock_version')
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
        default_vi_info = [ (x[1],x[3]) for x in vi_file_fields if x[0][0:3]=="VI_" ]
        for vi_key, vi_value in default_vi_info:
            self.cds_metadata.add([vi_value for i in range(nspec)], name=vi_key)


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
            with open(importlib.resources.files('prospect').joinpath("data", f"{line_category}_lines.txt")) as LINES:
                data = LINES.readlines()
            name, longname, wavelength, vacuum, major = zip(*[line.strip().split(',') for line in data if not line.startswith('#')])
            # wavelength = np.array([float(w) for w in wavelength])
            vacuum = [k == 'True' for k in vacuum]
            major = [k == 'True' for k in major]
            # vacuum_wavelength = wavelength.copy()
            vacuum_wavelength = np.array([float(w) if vacuum[i] else _airtovac(float(w)) for i, w in enumerate(wavelength)])
            line_data['restwave'].extend(vacuum_wavelength.tolist())
            line_data['plotwave'].extend((vacuum_wavelength * (1 + z)).tolist())
            line_data['name'].extend(name)
            line_data['longname'].extend(longname)
            line_data['plotname'].extend(name)
            line_data['emission'].extend([line_category == 'emission']*len(name))
            line_data['major'].extend(major)
            # line_array = np.genfromtxt(importlib.resources.files('prospect').joinpath("data", f"{line_category}_lines.txt"),
            #                            delimiter=",",
            #                            dtype=[("name", "|U20"),
            #                                   ("longname", "|U20"),
            #                                   ("wavelength", float),
            #                                   ("vacuum", bool),
            #                                   ("major", bool)],
            #                             encoding='utf-8')
            # vacuum_wavelengths = line_array['wavelength']
            # w, = np.where(line_array['vacuum']==False)
            # vacuum_wavelengths[w] = np.array([_airtovac(wave) for wave in line_array['wavelength'][w]])
            # line_data['restwave'].extend(vacuum_wavelengths)
            # line_data['plotwave'].extend(vacuum_wavelengths * (1+z))
            # line_data['name'].extend(line_array['name'])
            # line_data['longname'].extend(line_array['longname'])
            # line_data['plotname'].extend(line_array['name'])
            # emission_flag = True if line_category=='emission' else False
            # line_data['emission'].extend([emission_flag for row in line_array])
            # line_data['major'].extend(line_array['major'])

        self.cds_spectral_lines = ColumnDataSource(line_data)

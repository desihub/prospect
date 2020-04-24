# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
prospect.specutils
==================

Conversions of DESI & SDSS spectra into specutils-compatible objects.

"""
import os
import re

import numpy as np
import astropy.units as u
from astropy.nddata import InverseVariance
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from specutils import SpectrumList, Spectrum1D

from desiutil.depend import add_dependencies
from desiutil.io import encode_table

from desispec.maskbits import specmask
from desispec.resolution import Resolution
from desispec.io.util import fitsheader, native_endian, add_columns
from desispec.io.frame import read_frame
from desispec.io.fibermap import fibermap_comments
from .utils_specviewer import _coadd

class Spectra(SpectrumList):
    """Represents a grouping of spectra.

    This class contains an "extended" fibermap that has information about
    the night and exposure of each spectrum.  For each band, this class has
    the wavelength grid, flux, ivar, mask, and resolution arrays.

    Parameters
    ----------
    bands : :class:`list`
        List of strings used to identify the bands.
    wave : :class:`dict`
        Dictionary of arrays specifying the wavelength grid.
    flux : :class:`dict`
        Dictionary of arrays specifying the flux for each spectrum.
    ivar : :class:`dict`
        Dictionary of arrays specifying the inverse variance.
    mask : :class:`dict`, optional
        Dictionary of arrays specifying the bitmask.
    resolution_data : :class:`dict`, optional
        Dictionary of arrays specifying the block diagonal resolution matrix.
        The object for each band must be in one of the formats supported
        by the Resolution class constructor.
    fibermap
        Extended fibermap to use. If not specified, a fake one is created.
    meta : :class:`dict`, optional
        Dictionary of arbitrary properties.
    extra : :class:`dict`, optional
        Optional dictionary of dictionaries containing extra
        floating point arrays.  The top-level is a dictionary over bands
        and each value is a dictionary containing string keys and values
        which are arrays of the same size as the flux array.
    single : :class:`bool`, optional
        If ``True``, store data in memory as single precision.
    scores
        QA scores table.
    """
    def __init__(self, bands=[], wave={}, flux={}, ivar={}, mask=None, resolution_data=None,
        fibermap=None, meta=None, extra=None, single=False, scores=None):

        self._bands = tuple(bands)
        self._single = single
        self._ftype = np.float64
        if single:
            self._ftype = np.float32

        self._reset_properties()

        self.scores = scores

        if meta is None:
            self.meta = dict()
        elif isinstance(meta, fits.Header):
            self.meta = {'header': meta}
        else:
            self.meta = meta.copy()

        nspec = 0

        # check consistency of input dimensions
        for b in self._bands:
            if wave[b].ndim != 1:
                raise ValueError("Wavelength array for band {} should have shape (Nwave, ).".format(b))
            if flux[b].ndim != 2:
                raise ValueError("Flux array for band {} should have shape (Nspec, Nwave).".format(b))
            if flux[b].shape[1] != wave[b].shape[0]:
                raise ValueError("Flux array wavelength dimension for band {} does not match wavelength grid.".format(b))
            if nspec is None:
                nspec = flux[b].shape[0]
            if fibermap is not None:
                if len(fibermap) != flux[b].shape[0]:
                    raise ValueError("Flux array number of spectra for band {} does not match fibermap.".format(b))
            if ivar[b].shape != flux[b].shape:
                raise ValueError("Inverse variance array dimensions do not match flux for band {}.".format(b))
            if mask is not None:
                if mask[b].shape != flux[b].shape:
                    raise ValueError("Mask array dimensions do not match flux for band {}.".format(b))
                if mask[b].dtype not in (int, np.int64, np.int32, np.uint64, np.uint32):
                    raise ValueError("Bad mask type {}.".format(mask.dtype))
            if resolution_data is not None:
                if resolution_data[b].ndim != 3:
                    raise ValueError("Resolution array for band {} should have shape (Nspec, Ndiag, Nwave).".format(b))
                if resolution_data[b].shape[0] != flux[b].shape[0]:
                    raise ValueError("Resolution array spectrum dimension for band {} does not match flux.".format(b))
                if resolution_data[b].shape[2] != wave[b].shape[0]:
                    raise ValueError("Resolution array wavelength dimension for band {} does not match wavelength grid.".format(b))
            if extra is not None:
                for ex in extra[b].items():
                    if ex[1].shape != flux[b].shape:
                        raise ValueError("Extra arrays must have the same shape as the flux array.")

        if fibermap is not None:
            self.fibermap = fibermap.copy()
        else:
            self.fibermap = None

        # copy band-based data

        for b in self._bands:
            band_meta = dict()
            if mask is None:
                band_meta['mask'] = None
                bool_mask = None
            else:
                band_meta['mask'] = np.copy(mask[b])
                bool_mask = band_meta['mask'] != 0
            if resolution_data is None:
                band_meta['resolution_data'] = None
                band_meta['R'] = None
            else:
                band_meta['resolution_data'] = np.copy(resolution_data[b].astype(self._ftype))
                band_meta['R'] = np.array([Resolution(r) for r in resolution_data[b]])
            if extra is None:
                band_meta['extra'] = None
            else:
                band_meta['extra'] = dict()
                for k, v in extra[b].items():
                    band_meta['extra'][k] = np.copy(v.astype(self._ftype))
            self.append(Spectrum1D(spectral_axis=np.copy(wave[b].astype(self._ftype))*u.Angstrom,
                                   flux=np.copy(flux[b].astype(self._ftype))*u.Unit('10**-17 erg/(s cm2 Angstrom)'),
                                   uncertainty=InverseVariance(np.copy(ivar[b].astype(self._ftype))),
                                   mask=bool_mask,
                                   meta=band_meta))

    @property
    def bands(self):
        """
        (list): the list of valid bands.
        """
        return self._bands

    @property
    def ftype(self):
        """
        (numpy.dtype): the data type used for floating point numbers.
        """
        return self._ftype

    @property
    def wave(self):
        if self._wave is None:
            self._wave = dict()
            for band in self.bands:
                self._wave[band] = self.wavelength_grid(band)
        return self._wave

    @property
    def flux(self):
        if self._flux is None:
            self._flux = dict()
            for band in self.bands:
                self._flux[band] = self[self.bands.index(band)].flux
        return self._flux

    @property
    def ivar(self):
        if self._ivar is None:
            self._ivar = dict()
            for band in self.bands:
                self._ivar[band] = self[self.bands.index(band)].uncertainty
        return self._ivar

    @property
    def mask(self):
        if self._mask is None:
            self._mask = dict()
            for band in self.bands:
                self._mask[band] = self[self.bands.index(band)].meta['mask']
        return self._mask

    @property
    def resolution_data(self):
        if self._resolution_data is None:
            if self[0].meta['resolution_data'] is None:
                return None
            self._resolution_data = dict()
            for band in self.bands:
                self._resolution_data[band] = self[self.bands.index(band)].meta['resolution_data']
        return self._resolution_data

    @property
    def R(self):
        if self._R is None:
            if self[0].meta['R'] is None:
                return None
            self._R = dict()
            for band in self.bands:
                self._R[band] = self[self.bands.index(band)].meta['R']
        return self._R

    @property
    def extra(self):
        if self._extra is None:
            if self[0].meta['extra'] is None:
                return None
            self._extra = dict()
            for band in self.bands:
                self._extra[band] = self[self.bands.index(band)].meta['extra']
        return self._extra

    def _reset_properties(self):
        self._wave = None
        self._flux = None
        self._ivar = None
        self._mask = None
        self._resolution_data = None
        self._R = None
        self._extra = None

    def wavelength_grid(self, band):
        """
        Return the wavelength grid for a band.

        Args:
            band (str): the name of the band.

        Returns (array):
            an array containing the wavelength values.

        """
        if band not in self.bands:
            raise KeyError("{} is not a valid band.".format(band))
        return self.wave[band]

    def target_ids(self):
        """
        Return list of unique target IDs.

        The target IDs are sorted by the order that they first appear.

        Returns (array):
            an array of integer target IDs.
        """
        uniq, indices = np.unique(self.fibermap["TARGETID"], return_index=True)
        return uniq[indices.argsort()]


    def num_spectra(self):
        """
        Get the number of spectra contained in this group.

        Returns (int):
            Number of spectra contained in this group.
        """
        if self.fibermap is not None:
            return len(self.fibermap)
        else:
            return 0


    def num_targets(self):
        """
        Get the number of distinct targets.

        Returns (int):
            Number of unique targets with spectra in this object.
        """
        if self.fibermap is not None:
            return len(np.unique(self.fibermap["TARGETID"]))
        else:
            return 0


    def select(self, nights=None, bands=None, targets=None, fibers=None, invert=False):
        """
        Select a subset of the data.

        This filters the data based on a logical AND of the different
        criteria, optionally inverting that selection.

        Args:
            nights (list): optional list of nights to select.
            bands (list): optional list of bands to select.
            targets (list): optional list of target IDs to select.
            fibers (list): list/array of fiber indices to select.
            invert (bool): after combining all criteria, invert selection.

        Returns (Spectra):
            a new Spectra object containing the selected data.
        """
        keep_bands = None
        if bands is None:
            keep_bands = self.bands
        else:
            keep_bands = [ x for x in self.bands if x in bands ]
        if len(keep_bands) == 0:
            raise ValueError("No valid bands were selected!")

        keep_nights = None
        if nights is None:
            keep_nights = [ True for x in self.fibermap["NIGHT"] ]
        else:
            keep_nights = [ (x in nights) for x in self.fibermap["NIGHT"] ]
        if sum(keep_nights) == 0:
            raise ValueError("No valid nights were selected!")

        keep_targets = None
        if targets is None:
            keep_targets = [ True for x in self.fibermap["TARGETID"] ]
        else:
            keep_targets = [ (x in targets) for x in self.fibermap["TARGETID"] ]
        if sum(keep_targets) == 0:
            raise ValueError("No valid targets were selected!")

        keep_fibers = None
        if fibers is None:
            keep_fibers = [ True for x in self.fibermap["FIBER"] ]
        else:
            keep_fibers = [ (x in fibers) for x in self.fibermap["FIBER"] ]
        if sum(keep_fibers) == 0:
            raise ValueError("No valid fibers were selected!")

        keep_rows = [ (x and y and z) for x, y, z in zip(keep_nights, keep_targets, keep_fibers) ]
        if invert:
            keep_rows = [ not x for x in keep_rows ]

        keep = [ i for i, x in enumerate(keep_rows) if x ]
        if len(keep) == 0:
            raise ValueError("Selection has no spectra!")

        keep_wave = {}
        keep_flux = {}
        keep_ivar = {}
        keep_mask = None
        keep_res = None
        keep_extra = None
        if self.mask is not None:
            keep_mask = {}
        if self.resolution_data is not None:
            keep_res = {}
        if self.extra is not None:
            keep_extra = {}

        for b in keep_bands:
            keep_wave[b] = self.wave[b]
            keep_flux[b] = self.flux[b][keep,:]
            keep_ivar[b] = self.ivar[b][keep,:]
            if self.mask is not None:
                keep_mask[b] = self.mask[b][keep,:]
            if self.resolution_data is not None:
                keep_res[b] = self.resolution_data[b][keep,:,:]
            if self.extra is not None:
                keep_extra[b] = {}
                for k, v in self.extra[b].items():
                    keep_extra[b][k] = v[keep,:]

        return Spectra(keep_bands, keep_wave, keep_flux, keep_ivar,
                       mask=keep_mask, resolution_data=keep_res,
                       fibermap=self.fibermap[keep], meta=self.meta, extra=keep_extra,
                       single=self._single, scores=self.scores[keep])

    def update(self, other):
        """
        Overwrite or append new data.

        Given another Spectra object, compare the fibermap information with
        the existing one.  For spectra that already exist, overwrite existing
        data with the new values.  For spectra that do not exist, append that
        data to the end of the spectral data.

        Args:
            other (Spectra): the new data to add.

        Returns:
            nothing (object updated in place).

        """
        if not isinstance(other, Spectra):
            raise ValueError("New data has incorrect type!")

        # Does the other Spectra object have any data?

        if other.num_spectra() == 0:
            return

        # Do we have new bands to add?

        newbands = []
        for b in other.bands:
            if b not in self.bands:
                newbands.append(b)
            else:
                if not np.allclose(self.wave[b], other.wave[b]):
                    raise ValueError("Band {} has an incompatible wavelength grid.".format(b))

        bands = list(self.bands)
        bands.extend(newbands)

        # Are we adding mask data in this update?

        add_mask = False
        if other.mask is None:
            if self.mask is not None:
                raise ValueError("Existing spectra has a mask, cannot "
                                 "update it to a spectra with no mask.")
        else:
            if self.mask is None:
                add_mask = True

        # Are we adding resolution data in this update?

        ndiag = {}

        add_res = False
        if other.resolution_data is None:
            if self.resolution_data is not None:
                raise ValueError("Existing spectra has resolution data, cannot "
                                 "update it to a spectra with none.")
        else:
            if self.resolution_data is not None:
                for b in self.bands:
                    ndiag[b] = self.resolution_data[b].shape[1]
                for b in other.bands:
                    odiag = other.resolution_data[b].shape[1]
                    if b not in self.bands:
                        ndiag[b] = odiag
                    else:
                        if odiag != ndiag[b]:
                            raise ValueError("Resolution matrices for a"
                                             " given band must have the same dimensions.")
            else:
                add_res = True
                for b in other.bands:
                    ndiag[b] = other.resolution_data[b].shape[1]

        # Are we adding extra data in this update?

        add_extra = False
        if other.extra is None:
            if self.extra is not None:
                raise ValueError("Existing spectra has extra data, cannot "
                                 "update it to a spectra with none.")
        else:
            if self.extra is None:
                add_extra = True

        # Compute which targets / exposures are new

        nother = len(other.fibermap)
        exists = np.zeros(nother, dtype=np.int)

        indx_original = []

        if self.fibermap is not None:
            for r in range(nother):
                expid = other.fibermap[r]["EXPID"]
                fiber = other.fibermap[r]["FIBER"]
                for i, row in enumerate(self.fibermap):
                    if (expid == row["EXPID"]) and (fiber == row["FIBER"]):
                        indx_original.append(i)
                        exists[r] += 1

        if len(np.where(exists > 1)[0]) > 0:
            raise ValueError("Found duplicate spectra (same EXPID and FIBER) in the fibermap.")

        indx_exists = np.where(exists == 1)[0]
        indx_new = np.where(exists == 0)[0]

        # Make new data arrays of the correct size to hold both the old and
        # new data

        nupdate = len(indx_exists)
        nnew = len(indx_new)

        if self.fibermap is None:
            nold = 0
            newfmap = other.fibermap.copy()
        else:
            nold = len(self.fibermap)
            newfmap = encode_table(np.zeros( (nold + nnew, ),
                                   dtype=self.fibermap.dtype))

        if self.scores is None:
            if other.scores is None:
                newscores = None
            else:
                newscores = other.scores.copy()
        else:
            newscores = encode_table(np.zeros( (nold + nnew, ),
                                     dtype=self.scores.dtype))

        newwave = {}
        newflux = {}
        newivar = {}

        newmask = None
        if add_mask or self.mask is not None:
            newmask = {}

        newres = None
        newR = None
        if add_res or self.resolution_data is not None:
            newres = {}
            newR = {}

        newextra = None
        if add_extra or self.extra is not None:
            newextra = {}

        for b in bands:
            nwave = None
            if b in self.bands:
                nwave = self.wave[b].shape[0]
                newwave[b] = self.wave[b]
            else:
                nwave = other.wave[b].shape[0]
                newwave[b] = other.wave[b].astype(self._ftype)
            newflux[b] = np.zeros( (nold + nnew, nwave), dtype=self._ftype)
            newivar[b] = np.zeros( (nold + nnew, nwave), dtype=self._ftype)
            if newmask is not None:
                newmask[b] = np.zeros( (nold + nnew, nwave), dtype=np.uint32)
                newmask[b][:,:] = specmask["NODATA"]
            if newres is not None:
                newres[b] = np.zeros( (nold + nnew, ndiag[b], nwave), dtype=self._ftype)
            if newextra is not None:
                newextra[b] = {}

        # Copy the old data

        if nold > 0:
            # We have some data (i.e. we are not starting with an empty Spectra)
            newfmap[:nold] = self.fibermap
            if newscores is not None:
                newscores[:nold] = self.scores

            for b in self.bands:
                newflux[b][:nold,:] = self.flux[b]
                newivar[b][:nold,:] = self.ivar[b]
                if self.mask is not None:
                    newmask[b][:nold,:] = self.mask[b]
                elif add_mask:
                    newmask[b][:nold,:] = 0
                if self.resolution_data is not None:
                    newres[b][:nold,:,:] = self.resolution_data[b]
                if self.extra is not None:
                    for ex in self.extra[b].items():
                        newextra[b][ex[0]] = np.zeros( newflux[b].shape,
                            dtype=self._ftype)
                        newextra[b][ex[0]][:nold,:] = ex[1]

        # Update existing spectra

        for i, s in enumerate(indx_exists):
            row = indx_original[i]
            for b in other.bands:
                newflux[b][row,:] = other.flux[b][s,:].astype(self._ftype)
                newivar[b][row,:] = other.ivar[b][s,:].astype(self._ftype)
                if other.mask is not None:
                    newmask[b][row,:] = other.mask[b][s,:]
                else:
                    newmask[b][row,:] = 0
                if other.resolution_data is not None:
                    newres[b][row,:,:] = other.resolution_data[b][s,:,:].astype(self._ftype)
                if other.extra is not None:
                    for ex in other.extra[b].items():
                        if ex[0] not in newextra[b]:
                            newextra[b][ex[0]] = np.zeros(newflux[b].shape,
                                dtype=self._ftype)
                        newextra[b][ex[0]][row,:] = ex[1][s,:].astype(self._ftype)

        # Append new spectra

        if nnew > 0:
            newfmap[nold:] = other.fibermap[indx_new]
            if newscores is not None:
                newscores[nold:] = other.scores[indx_new]

            for b in other.bands:
                newflux[b][nold:,:] = other.flux[b][indx_new].astype(self._ftype)
                newivar[b][nold:,:] = other.ivar[b][indx_new].astype(self._ftype)
                if other.mask is not None:
                    newmask[b][nold:,:] = other.mask[b][indx_new]
                else:
                    newmask[b][nold:,:] = 0
                if other.resolution_data is not None:
                    newres[b][nold:,:,:] = other.resolution_data[b][indx_new].astype(self._ftype)
                if other.extra is not None:
                    for ex in other.extra[b].items():
                        if ex[0] not in newextra[b]:
                            newextra[b][ex[0]] = np.zeros(newflux[b].shape,
                                dtype=self._ftype)
                        newextra[b][ex[0]][nold:,:] = ex[1][indx_new].astype(self._ftype)

        # Swap data into place

        self._bands = bands
        self.fibermap = newfmap
        self.scores = newscores
        self._reset_properties()
        for i, b in enumerate(self._bands):
            band_meta = dict()
            if newmask is None:
                band_meta['mask'] = None
                bool_mask = None
            else:
                band_meta['mask'] = newmask[b]
                bool_mask = band_meta['mask'] != 0
            if newres is None:
                band_meta['resolution_data'] = None
                band_meta['R'] = None
            else:
                band_meta['resolution_data'] = newres[b]
                band_meta['R'] = np.array([Resolution(r) for r in newres[b]])
            if newextra is None:
                band_meta['extra'] = None
            else:
                band_meta['extra'] = dict()
                for k, v in newextra[b].items():
                    band_meta['extra'][k] = v
            s = Spectrum1D(spectral_axis=newwave[b]*u.Angstrom,
                           flux=newflux[b]*u.Unit('10**-17 erg/(s cm2 Angstrom)'),
                           uncertainty=InverseVariance(newivar),
                           mask=bool_mask,
                           meta=band_meta)
            try:
                self[i] = s
            except IndexError:
                self.append(s)
        return


def write_spectra(outfile, spec, units=None):
    """
    Write Spectra object to FITS file.

    This places the metadata into the header of the (empty) primary HDU.
    The first extension contains the fibermap, and then HDUs are created for
    the different data arrays for each band.

    Floating point data is converted to 32 bits before writing.

    Args:
        outfile (str): path to write
        spec (Spectra): the object containing the data
        units (str): optional string to use for the BUNIT key of the flux
            HDUs for each band.

    Returns:
        The absolute path to the file that was written.

    """

    outfile = os.path.abspath(outfile)

    # Create the parent directory, if necessary.
    dir, base = os.path.split(outfile)
    if not os.path.exists(dir):
        os.makedirs(dir)

    # Create HDUs from the data
    all_hdus = fits.HDUList()

    # metadata goes in empty primary HDU
    hdr = fitsheader(spec.meta)
    add_dependencies(hdr)

    all_hdus.append(fits.PrimaryHDU(header=hdr))

    # Next is the fibermap
    fmap = spec.fibermap.copy()
    fmap.meta["EXTNAME"] = "FIBERMAP"
    hdu = fits.convenience.table_to_hdu(fmap)

    # Add comments for fibermap columns.
    for i, colname in enumerate(fmap.dtype.names):
        if colname in fibermap_comments:
            key = "TTYPE{}".format(i+1)
            name = hdu.header[key]
            assert name == colname
            comment = fibermap_comments[name]
            hdu.header[key] = (name, comment)
        else:
            pass
            #print('Unknown comment for {}'.format(colname))

    all_hdus.append(hdu)

    # Now append the data for all bands

    for band in spec.bands:
        hdu = fits.ImageHDU(name="{}_WAVELENGTH".format(band.upper()))
        hdu.header["BUNIT"] = "Angstrom"
        hdu.data = spec.wave[band].astype("f8")
        all_hdus.append(hdu)

        hdu = fits.ImageHDU(name="{}_FLUX".format(band.upper()))
        if units is None:
            hdu.header["BUNIT"] = "10**-17 erg/(s cm2 Angstrom)"
        else:
            hdu.header["BUNIT"] = units
        hdu.data = spec.flux[band].astype("f4")
        all_hdus.append(hdu)

        hdu = fits.ImageHDU(name="{}_IVAR".format(band.upper()))
        if units is None:
            hdu.header["BUNIT"] = '10**+34 (s2 cm4 Angstrom2) / erg2'
        else:
            hdu.header["BUNIT"] = ((u.Unit(units, format='fits'))**-2).to_string('fits')
        hdu.data = spec.ivar[band].astype("f4")
        all_hdus.append(hdu)

        if spec.mask is not None:
            # hdu = fits.CompImageHDU(name="{}_MASK".format(band.upper()))
            hdu = fits.ImageHDU(name="{}_MASK".format(band.upper()))
            hdu.data = spec.mask[band].astype(np.uint32)
            all_hdus.append(hdu)

        if spec.resolution_data is not None:
            hdu = fits.ImageHDU(name="{}_RESOLUTION".format(band.upper()))
            hdu.data = spec.resolution_data[band].astype("f4")
            all_hdus.append(hdu)

        if spec.extra is not None:
            for ex in spec.extra[band].items():
                hdu = fits.ImageHDU(name="{}_{}".format(band.upper(), ex[0]))
                hdu.data = ex[1].astype("f4")
                all_hdus.append(hdu)

    if spec.scores is not None :
        scores_tbl = encode_table(spec.scores)  #- unicode -> bytes
        scores_tbl.meta['EXTNAME'] = 'SCORES'
        all_hdus.append( fits.convenience.table_to_hdu(scores_tbl) )
        if spec.scores_comments is not None : # add comments in header
            hdu=all_hdus['SCORES']
            for i in range(1,999):
                key = 'TTYPE'+str(i)
                if key in hdu.header:
                    value = hdu.header[key]
                    if value in spec.scores_comments.keys() :
                        hdu.header[key] = (value, spec.scores_comments[value])

    all_hdus.writeto("{}.tmp".format(outfile), overwrite=True, checksum=True)
    os.rename("{}.tmp".format(outfile), outfile)

    return outfile


def read_spectra(infile, single=False, coadd=False):
    """Read Spectra object from FITS file.

    This reads data written by the write_spectra function.  A new Spectra
    object is instantiated and returned.

    Args:
        infile (str): path to read
        single (bool): if True, keep spectra as single precision in memory.
        coadd (bool): if True, coadd all spectra from the same targetid.

    Returns (Spectra):
        The object containing the data read from disk.

    """

    ftype = np.float64
    if single:
        ftype = np.float32

    infile = os.path.abspath(infile)
    if not os.path.isfile(infile):
        raise FileNotFoundError("{} is not a file!".format(infile))

    # initialize data objects
    bands = []
    fmap = None
    wave = None
    flux = None
    ivar = None
    mask = None
    res = None
    extra = None
    scores = None

    with fits.open(infile, mode="readonly") as hdulist:
        nhdu = len(hdulist)

        # load the metadata.
        meta = hdulist[0].header

        # For efficiency, go through the HDUs in disk-order.  Use the
        # extension name to determine where to put the data.  We don't
        # explicitly copy the data, since that will be done when constructing
        # the Spectra object.

        for h in range(1, nhdu):
            name = hdulist[h].header["EXTNAME"]
            if name == "FIBERMAP":
                fmap = encode_table(Table(hdulist[h].data, copy=True).as_array())
            elif name == "SCORES":
                scores = encode_table(Table(hdulist[h].data, copy=True).as_array())
            else:
                # Find the band based on the name
                mat = re.match(r"(.*)_(.*)", name)
                if mat is None:
                    raise RuntimeError("FITS extension name {} does not contain the band".format(name))
                band = mat.group(1).lower()
                type = mat.group(2)
                if band not in bands:
                    bands.append(band)
                if type == "WAVELENGTH":
                    if wave is None:
                        wave = {}
                    wave[band] = native_endian(hdulist[h].data.astype(ftype))
                elif type == "FLUX":
                    if flux is None:
                        flux = {}
                    flux[band] = native_endian(hdulist[h].data.astype(ftype))
                elif type == "IVAR":
                    if ivar is None:
                        ivar = {}
                    ivar[band] = native_endian(hdulist[h].data.astype(ftype))
                elif type == "MASK":
                    if mask is None:
                        mask = {}
                    mask[band] = native_endian(hdulist[h].data.astype(np.uint32))
                elif type == "RESOLUTION":
                    if res is None:
                        res = {}
                    res[band] = native_endian(hdulist[h].data.astype(ftype))
                else:
                    # this must be an "extra" HDU
                    if extra is None:
                        extra = {}
                    if band not in extra:
                        extra[band] = {}
                    extra[band][type] = native_endian(hdulist[h].data.astype(ftype))

    if coadd:
        uniq, indices = np.unique(fmap["TARGETID"], return_index=True)
        targetids = uniq[indices.argsort()]
        ntargets = len(targetids)
        cwave = dict()
        cflux = dict()
        civar = dict()
        crdat = dict()
        cmask = dict()
        for channel in bands:
            cwave[channel] = wave[channel].copy()
            nwave = len(cwave[channel])
            cflux[channel] = np.zeros((ntargets, nwave))
            civar[channel] = np.zeros((ntargets, nwave))
            ndiag = res[channel].shape[1]
            crdat[channel] = np.zeros((ntargets, ndiag, nwave))
            cmask[channel] = np.zeros((ntargets, nwave), dtype=mask[channel].dtype)
        #- Loop over targets, coadding all spectra for each target
        fibermap = Table(dtype=fmap.dtype)
        for i, targetid in enumerate(targetids):
            ii = np.where(fmap['TARGETID'] == targetid)[0]
            fibermap.add_row(fmap[ii[0]])
            for channel in bands:
                if len(ii) > 1:
                    outwave, outflux, outivar, outrdat = _coadd(
                        wave[channel],
                        flux[channel][ii],
                        ivar[channel][ii],
                        res[channel][ii]
                        )
                    outmask = mask[channel][ii[0]]
                    for j in range(1, len(ii)):
                        outmask |= mask[channel][ii[j]]
                else:
                    outwave, outflux, outivar, outrdat = (
                        wave[channel],
                        flux[channel][ii[0]],
                        ivar[channel][ii[0]],
                        res[channel][ii[0]]
                        )
                    outmask = mask[channel][ii[0]]

                cflux[channel][i] = outflux
                civar[channel][i] = outivar
                crdat[channel][i] = outrdat
                cmask[channel][i] = outmask

        return Spectra(bands, cwave, cflux, civar, mask=cmask, resolution_data=crdat,
                       fibermap=fibermap, meta=meta, extra=extra, single=single,
                       scores=scores)

    # Construct the Spectra object from the data.  If there are any
    # inconsistencies in the sizes of the arrays read from the file,
    # they will be caught by the constructor.

    return Spectra(bands, wave, flux, ivar, mask=mask, resolution_data=res,
                   fibermap=fmap, meta=meta, extra=extra, single=single,
                   scores=scores)


def read_spPlate(filename):
    """Read a SDSS spPlate file.

    Parameters
    ----------
    filename : :class:`str`
        Name of the spPlate file.

    Returns
    -------
    Spectrum1D
        The spectra.
    """
    with fits.open(filename) as hdulist:
        header = hdulist[0].header
        meta = {'header': header}
        try:
            flux_unit = u.Unit(hdulist[0].header['BUNIT'])
        except ValueError:
            flux_unit = u.Unit('1e-17 erg / (Angstrom cm2 s)')
        flux = hdulist[0].data * flux_unit
        wcs = WCS(header)
        dispersion_unit = u.Unit('Angstrom')
        dispersion = 10**wcs.all_pix2world(np.vstack((np.arange(flux.shape[1]),
                                                      np.zeros((flux.shape[1],)))).T,
                                           0)[:, 0]
        try:
            uncertainty_unit = u.Unit(hdulist[1].header['BUNIT'])
        except ValueError:
            uncertainty_unit = u.Unit('1e+34 (Angstrom2 cm4 s2) / erg2')
        uncertainty = InverseVariance(hdulist[1].data * uncertainty_unit)
        mask = hdulist[2].data != 0
        meta['plugmap'] = Table.read(hdulist[5])
    return Spectrum1D(flux=flux, spectral_axis=dispersion*dispersion_unit,
                      uncertainty=uncertainty, meta=meta, mask=mask)


def read_spZbest(filename):
    """Read a SDSS spZbest file.

    Parameters
    ----------
    filename : :class:`str`
        Name of the spZbest file.

    Returns
    -------
    tuple
        A Table containing the redshift values and a Spectrum1D object containing
        the best-fit models.
    """
    with fits.open(filename) as hdulist:
        header = hdulist[0].header
        meta = {'header': header}
        redshifts = Table.read(hdulist[1])
        flux_unit = u.Unit('1e-17 erg / (Angstrom cm2 s)')
        flux = hdulist[2].data * flux_unit
        dispersion = 10**(header['CRVAL1'] +
                          header['CD1_1'] * np.arange(hdulist[2].header['NAXIS1'],
                                                      dtype=hdulist[2].data.dtype))
        dispersion_unit = u.Unit('Angstrom')
        models = Spectrum1D(flux=flux, spectral_axis=dispersion*dispersion_unit,
                            meta=meta)
    return redshifts, models


def read_frame_as_spectra(filename, night=None, expid=None, band=None, single=False):
    """
    Read a FITS file containing a Frame and return a Spectra.

    A Frame file is very close to a Spectra object (by design), and
    only differs by missing the NIGHT and EXPID in the fibermap, as
    well as containing only one band of data.

    Args:
        infile (str): path to read

    Options:
        night (int): the night value to use for all rows of the fibermap.
        expid (int): the expid value to use for all rows of the fibermap.
        band (str): the name of this band.
        single (bool): if True, keep spectra as single precision in memory.

    Returns (Spectra):
        The object containing the data read from disk.

    """
    fr = read_frame(filename)
    if fr.fibermap is None:
        raise RuntimeError("reading Frame files into Spectra only supported if a fibermap exists")

    nspec = len(fr.fibermap)

    if band is None:
        band = fr.meta['camera'][0]

    if night is None:
        night = fr.meta['night']

    if expid is None:
        expid = fr.meta['expid']

    fmap = np.asarray(fr.fibermap.copy())
    fmap = add_columns(fmap,
                       ['NIGHT', 'EXPID', 'TILEID'],
                       [np.int32(night), np.int32(expid), np.int32(fr.meta['TILEID'])],
                       )

    fmap = encode_table(fmap)

    bands = [ band ]

    mask = None
    if fr.mask is not None:
        mask = {band : fr.mask}

    res = None
    if fr.resolution_data is not None:
        res = {band : fr.resolution_data}

    extra = None
    if fr.chi2pix is not None:
        extra = {band : {"CHI2PIX" : fr.chi2pix}}

    spec = Spectra(bands, {band : fr.wave}, {band : fr.flux}, {band : fr.ivar},
        mask=mask, resolution_data=res, fibermap=fmap, meta=fr.meta,
        extra=extra, single=single, scores=fr.scores)

    return spec

# EA - TEMPORARY modification to desi spectra.update() function
# Does NOT rely on fibermap information to merge data : always append new data
# Reason : EXPID is not available in coadd spectra

import numpy as np
import desispec.spectra
from desiutil.io import encode_table
from desispec.maskbits import specmask
from desispec.resolution import Resolution

def myspecupdate(spectra_in, other) :
    """
    Args:
        other (Spectra): the new data to add.
    Returns:
        nothing (object "spectra" updated in place).
    """

    # Does the other Spectra object have any data?

    if other.num_spectra() == 0:
        return

    # Do we have new bands to add?

    newbands = []
    for b in other.bands:
        if b not in spectra_in.bands:
            newbands.append(b)
        else:
            if not np.allclose(spectra_in.wave[b], other.wave[b]):
                raise RuntimeError("band {} has an incompatible wavelength grid".format(b))

    bands = list(spectra_in.bands)
    bands.extend(newbands)

    # Are we adding mask data in this update?

    add_mask = False
    if other.mask is None:
        if spectra_in.mask is not None:
            raise RuntimeError("existing spectra has a mask, cannot "
                "update it to a spectra with no mask")
    else:
        if spectra_in.mask is None:
            add_mask = True

    # Are we adding resolution data in this update?

    ndiag = {}

    add_res = False
    if other.resolution_data is None:
        if spectra_in.resolution_data is not None:
            raise RuntimeError("existing spectra has resolution data, cannot "
                "update it to a spectra with none")
    else:
        if spectra_in.resolution_data is not None:
            for b in spectra_in.bands:
                ndiag[b] = spectra_in.resolution_data[b].shape[1]
            for b in other.bands:
                odiag = other.resolution_data[b].shape[1]
                if b not in spectra_in.bands:
                    ndiag[b] = odiag
                else:
                    if odiag != ndiag[b]:
                        raise RuntimeError("Resolution matrices for a"
                            " given band must have the same dimensoins")
        else:
            add_res = True
            for b in other.bands:
                ndiag[b] = other.resolution_data[b].shape[1]

    # Are we adding extra data in this update?

    add_extra = False
    if other.extra is None:
        if spectra_in.extra is not None:
            raise RuntimeError("existing spectra has extra data, cannot "
                "update it to a spectra with none")
    else:
        if spectra_in.extra is None:
            add_extra = True

    # Compute which targets / exposures are new

    nother = len(other.fibermap)
    exists = np.zeros(nother, dtype=np.int)

    indx_original = []

    # EA modif :
    check_exists = True
    if ( (spectra_in.fibermap is None)
        or ("EXPID" not in spectra_in.fibermap.keys())
        or ("EXPID" not in other.fibermap.keys())
        or ("FIBER" not in spectra_in.fibermap.keys())
        or ("FIBER" not in other.fibermap.keys()) ) :
        check_exists = False
    if check_exists : 
        for r in range(nother):
            expid = other.fibermap[r]["EXPID"]
            fiber = other.fibermap[r]["FIBER"]
            for i, row in enumerate(spectra_in.fibermap):
                if (expid == row["EXPID"]) and (fiber == row["FIBER"]):
                    indx_original.append(i)
                    exists[r] += 1

    if len(np.where(exists > 1)[0]) > 0:
        raise RuntimeError("found duplicate spectra (same EXPID and FIBER) in the fibermap")

    indx_exists = np.where(exists == 1)[0]
    indx_new = np.where(exists == 0)[0]

    # Make new data arrays of the correct size to hold both the old and 
    # new data

    nupdate = len(indx_exists)
    nnew = len(indx_new)

    if spectra_in.fibermap is None:
        nold = 0
        newfmap = other.fibermap.copy()
    else:
        nold = len(spectra_in.fibermap)
        newfmap = encode_table(np.zeros( (nold + nnew, ),
                               dtype=spectra_in.fibermap.dtype))

    newwave = {}
    newflux = {}
    newivar = {}

    newmask = None
    if add_mask or spectra_in.mask is not None:
        newmask = {}

    newres = None
    newR = None
    if add_res or spectra_in.resolution_data is not None:
        newres = {}
        newR = {}

    newextra = None
    if add_extra or spectra_in.extra is not None:
        newextra = {}

    for b in bands:
        nwave = None
        if b in spectra_in.bands:
            nwave = spectra_in.wave[b].shape[0]
            newwave[b] = spectra_in.wave[b]
        else:
            nwave = other.wave[b].shape[0]
            newwave[b] = other.wave[b].astype(spectra_in._ftype)
        newflux[b] = np.zeros( (nold + nnew, nwave), dtype=spectra_in._ftype)
        newivar[b] = np.zeros( (nold + nnew, nwave), dtype=spectra_in._ftype)
        if newmask is not None:
            newmask[b] = np.zeros( (nold + nnew, nwave), dtype=np.uint32)
            newmask[b][:,:] = specmask["NODATA"]
        if newres is not None:
            newres[b] = np.zeros( (nold + nnew, ndiag[b], nwave), dtype=spectra_in._ftype)
        if newextra is not None:
            newextra[b] = {}

    # Copy the old data

    if nold > 0:
        # We have some data (i.e. we are not starting with an empty Spectra)
        newfmap[:nold] = spectra_in.fibermap

        for b in spectra_in.bands:
            newflux[b][:nold,:] = spectra_in.flux[b]
            newivar[b][:nold,:] = spectra_in.ivar[b]
            if spectra_in.mask is not None:
                newmask[b][:nold,:] = spectra_in.mask[b]
            elif add_mask:
                newmask[b][:nold,:] = 0
            if spectra_in.resolution_data is not None:
                newres[b][:nold,:,:] = spectra_in.resolution_data[b]
            if spectra_in.extra is not None:
                for ex in spectra_in.extra[b].items():
                    newextra[b][ex[0]] = np.zeros( newflux[b].shape,
                        dtype=spectra_in._ftype)
                    newextra[b][ex[0]][:nold,:] = ex[1]

    # Update existing spectra

    for i, s in enumerate(indx_exists):
        row = indx_original[i]
        for b in other.bands:
            newflux[b][row,:] = other.flux[b][s,:].astype(spectra_in._ftype)
            newivar[b][row,:] = other.ivar[b][s,:].astype(spectra_in._ftype)
            if other.mask is not None:
                newmask[b][row,:] = other.mask[b][s,:]
            else:
                newmask[b][row,:] = 0
            if other.resolution_data is not None:
                newres[b][row,:,:] = other.resolution_data[b][s,:,:].astype(spectra_in._ftype)
            if other.extra is not None:
                for ex in other.extra[b].items():
                    if ex[0] not in newextra[b]:
                        newextra[b][ex[0]] = np.zeros(newflux[b].shape,
                            dtype=spectra_in._ftype)
                    newextra[b][ex[0]][row,:] = ex[1][s,:].astype(spectra_in._ftype)

    # Append new spectra

    if nnew > 0:
        newfmap[nold:] = other.fibermap[indx_new]

        for b in other.bands:
            newflux[b][nold:,:] = other.flux[b][indx_new].astype(spectra_in._ftype)
            newivar[b][nold:,:] = other.ivar[b][indx_new].astype(spectra_in._ftype)
            if other.mask is not None:
                newmask[b][nold:,:] = other.mask[b][indx_new]
            else:
                newmask[b][nold:,:] = 0
            if other.resolution_data is not None:
                newres[b][nold:,:,:] = other.resolution_data[b][indx_new].astype(spectra_in._ftype)
            if other.extra is not None:
                for ex in other.extra[b].items():
                    if ex[0] not in newextra[b]:
                        newextra[b][ex[0]] = np.zeros(newflux[b].shape,
                            dtype=spectra_in._ftype)
                    newextra[b][ex[0]][nold:,:] = ex[1][indx_new].astype(spectra_in._ftype)

    # Update all sparse resolution matrices

    for b in bands:
        if newres is not None:
            newR[b] = np.array( [ Resolution(r) for r in newres[b] ] )

    # Swap data into place

    spectra_in._bands = bands
    spectra_in.wave = newwave
    spectra_in.fibermap = newfmap
    spectra_in.flux = newflux
    spectra_in.ivar = newivar
    spectra_in.mask = newmask
    spectra_in.resolution_data = newres
    spectra_in.R = newR
    spectra_in.extra = newextra

    return spectra_in

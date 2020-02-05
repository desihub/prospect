
# EA - June 2019. TEMPORARY modification to desi spectra.select() function
# to include expid-based selection + indices-based selection
# changes : fct name+args ; self=>thespec ; expid/indices-based selection + final selection

import desispec.spectra

def myspecselect(thespec, nights=None, bands=None, targets=None, fibers=None, expids=None, indices=None, invert=False):
    """
    Select a subset of the data.
    This filters the data based on a logical AND of the different
    criteria, optionally inverting that selection.
    Args:
        nights (list): optional list of nights to select.
        bands (list): optional list of bands to select.
        targets (list): optional list of target IDs to select.
        fibers (list): list/array of fiber indices to select.
        ADDED=> expids (list): list/array of individual exposures to select.      
        ADDED =>indices (list) : list of raw (arbitrary) indices in the Spectra object to select. 
        invert (bool): after combining all criteria, invert selection.
    Returns (Spectra):
        a new Spectra object containing the selected data.
    """
    keep_bands = None
    if bands is None:
        keep_bands = thespec.bands
    else:
        keep_bands = [ x for x in thespec.bands if x in bands ]
    if len(keep_bands) == 0:
        raise RuntimeError("no valid bands were selected!")

    keep_nights = None
    if nights is None:
        keep_nights = [ True for x in range(thespec.num_spectra()) ]
    else:
        keep_nights = [ (x in nights) for x in thespec.fibermap["NIGHT"] ]
    if sum(keep_nights) == 0:
        raise RuntimeError("no valid nights were selected!")

    keep_targets = None
    if targets is None:
        keep_targets = [ True for x in range(thespec.num_spectra()) ]
    else:
        keep_targets = [ (x in targets) for x in thespec.fibermap["TARGETID"] ]
    if sum(keep_targets) == 0:
        raise RuntimeError("no valid targets were selected!")

    keep_fibers = None
    if fibers is None:
        keep_fibers = [ True for x in range(thespec.num_spectra()) ]
    else:
        keep_fibers = [ (x in fibers) for x in thespec.fibermap["FIBER"] ]
    if sum(keep_fibers) == 0:
        raise RuntimeError("no valid fibers were selected!")

    keep_expids = None
    if expids is None:
        keep_expids = [ True for x in range(thespec.num_spectra()) ]
    else:
        keep_expids = [ (x in expids) for x in thespec.fibermap["EXPID"] ]
    if sum(keep_expids) == 0:
        raise RuntimeError("no valid expids were selected!")

    keep_indices = None
    if indices is None:
        keep_indices = [ True for x in range(thespec.num_spectra()) ]
    else:
        keep_indices = [ (x in indices) for x in range(thespec.num_spectra()) ]
    if sum(keep_indices) == 0:
        raise RuntimeError("no valid indices were selected!")

    keep_rows = [ (x and y and z and t and u) for x, y, z, t, u in zip(keep_nights, keep_targets, keep_fibers, keep_expids, keep_indices) ]
    if invert:
        keep_rows = [ not x for x in keep_rows ]

    keep = [ i for i, x in enumerate(keep_rows) if x ]
    if len(keep) == 0:
        raise RuntimeError("selection has no spectra")

    keep_wave = {}
    keep_flux = {}
    keep_ivar = {}
    keep_mask = None
    keep_res = None
    keep_extra = None
    if thespec.mask is not None:
        keep_mask = {}
    if thespec.resolution_data is not None:
        keep_res = {}
    if thespec.extra is not None:
        keep_extra = {}

    for b in keep_bands:
        keep_wave[b] = thespec.wave[b]
        keep_flux[b] = thespec.flux[b][keep,:]
        keep_ivar[b] = thespec.ivar[b][keep,:]
        if thespec.mask is not None:
            keep_mask[b] = thespec.mask[b][keep,:]
        if thespec.resolution_data is not None:
            keep_res[b] = thespec.resolution_data[b][keep,:,:]
        if thespec.extra is not None:
            keep_extra[b] = {}
            for ex in thespec.extra[b].items():
                keep_extra[b][ex[0]] = ex[1][keep,:]

    ret = desispec.spectra.Spectra(keep_bands, keep_wave, keep_flux, keep_ivar, 
        mask=keep_mask, resolution_data=keep_res, 
        fibermap=thespec.fibermap[keep], meta=thespec.meta, extra=keep_extra,
        single=thespec._single, scores=thespec.scores[keep])

    return ret

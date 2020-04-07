// Utility javascript function used in CustomJS callbacks
//  (Standalone module)

function adapt_plotrange(pmin, pmax, data) {
    // Used to adapt plotting range to a data vector
    
    // copy before sorting to not impact original, and filter out NaN
    var dx = data.slice().filter(Boolean)
    dx.sort()
    var imin = Math.floor(pmin * dx.length)
    var imax = Math.floor(pmax * dx.length)
    
    return [dx[imin], dx[imax]]
}


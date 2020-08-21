// Utility javascript function used in CustomJS callbacks
//  (Standalone module)

function adapt_plotrange(pmin, pmax, data) {
    // Used to adapt plotting range to a data vector
    // copy before sorting to not impact original, and filter out NaN
    var dx = data.slice().filter(isFinite);
    dx.sort(function(a, b){return a - b});
    var imin = Math.floor(pmin * dx.length);
    var imax = Math.floor(pmax * dx.length);
    // console.log("imin = " + imin + "; imax = " + imax);
    // console.log("dx[imin] = " + dx[imin] + "; dx[imax] = " + dx[imax]);
    // console.log("dx[0] = " + dx[0] + "; dx[dx.length-1] = " + dx[dx.length-1]);
    return [dx[imin], dx[imax]];
}

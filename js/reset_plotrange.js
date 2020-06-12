// CustomJS, callback for the "Reset X-Y range" button
//  - Requires to include function in: adapt_plotrange.js
//  - args = fig, xmin, xmax, spectra

// x-range : use fixed x-range determined once for all
fig.x_range.start = xmin
fig.x_range.end = xmax

var ymin = 0.0
var ymax = 0.0
for (var i=0; i<spectra.length; i++) {
    var data = spectra[i].data
    var tmp = adapt_plotrange(0.01, 0.99, data['plotflux'])
    ymin = Math.min(ymin, tmp[0])
    ymax = Math.max(ymax, tmp[1])
}
if(ymin<0) {
    fig.y_range.start = ymin * 1.4
} else {
    fig.y_range.start = ymin * 0.6
}
fig.y_range.end = ymax * 1.4

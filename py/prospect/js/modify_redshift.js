// CustomJS, callback for the "Redshift value" button
//  - args = spectra, coaddcam_spec, model, othermodel, metadata, ispectrumslider,
//           zslider, dzslider, z_input, widgetinfos, line_data, lines,
//           line_labels, zlines, zline_labels, overlap_waves, overlap_bands, fig.

var z = parseFloat(z_input.value)
if ( z >=-0.1 && z <= 5.0 ) {
    // update zsliders only if needed (avoid recursive call)
    z_input.value = parseFloat(z_input.value).toFixed(4)
    var z1 = Math.floor(z*100) / 100
    var z2 = z-z1
    // use a number slightly smaller than 0.01 to avoid rounding issue
    if ( Math.abs(z1-zslider.value) >= 0.00999999) zslider.value = parseFloat(parseFloat(z1).toFixed(2))
    if ( Math.abs(z2-dzslider.value) >= 0.00009999) dzslider.value = parseFloat(parseFloat(z2).toFixed(4))
} else {
    if (z_input.value < -0.1) z_input.value = (-0.1).toFixed(4)
    if (z_input.value > 5) z_input.value = (5.0).toFixed(4)
}
z = parseFloat(z_input.value)

var line_restwave = line_data.data['restwave']
var waveshift_lines = (widgetinfos.data['waveframe_active'][0] == 0) ? 1+z : 1 ;
var waveshift_spec = (widgetinfos.data['waveframe_active'][0] == 0) ? 1 : 1/(1+z) ;

for(var i=0; i<line_restwave.length; i++) {
    lines[i].location = line_restwave[i] * waveshift_lines
    line_labels[i].x = line_restwave[i] * waveshift_lines
    zlines[i].location = line_restwave[i] * waveshift_lines
    zline_labels[i].x = line_restwave[i] * waveshift_lines
}
if (overlap_bands.length>0) {
    for (var i=0; i<overlap_bands.length; i++) {
        overlap_bands[i].left = overlap_waves[i][0] * waveshift_spec
        overlap_bands[i].right = overlap_waves[i][1] * waveshift_spec
    }
}

function shift_plotwave(cds_spec, waveshift) {
    var data = cds_spec.data
    var origwave = data['origwave']
    var plotwave = data['plotwave']
    if ( plotwave[0] != origwave[0] * waveshift ) { // Avoid redo calculation if not needed
        for (var j=0; j<plotwave.length; j++) {
            plotwave[j] = origwave[j] * waveshift ;
        }
        cds_spec.change.emit()
    }
}

for(var i=0; i<spectra.length; i++) {
    shift_plotwave(spectra[i], waveshift_spec)
}
if (coaddcam_spec) shift_plotwave(coaddcam_spec, waveshift_spec)

// Update model wavelength array
// NEW : don't shift model if othermodel is there
if (othermodel) {
    var zref = othermodel.data['zref'][0]
    var waveshift_model = (widgetinfos.data['waveframe_active'][0] == 0) ? (1+z)/(1+zref) : 1/(1+zref) ;
    shift_plotwave(othermodel, waveshift_model)
} else if (model) {
    var zfit = 0.0
    if(metadata.data['Z'] !== undefined) {
        zfit = metadata.data['Z'][ispectrumslider.value]
    }
    var waveshift_model = (widgetinfos.data['waveframe_active'][0] == 0) ? (1+z)/(1+zfit) : 1/(1+zfit) ;
    shift_plotwave(model, waveshift_model)
}


// CustomJS, callback for the "Redshift value" button
//  - Requires to include functions in shift_wave.js
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
    if ( Math.abs(z1-zslider.value) >= 0.00999999 ) zslider.value = parseFloat(parseFloat(z1).toFixed(2))
    if ( Math.abs(z2-dzslider.value) >= 0.00009999 ) dzslider.value = parseFloat(parseFloat(z2).toFixed(4))
} else {
    if (z_input.value < -0.1) z_input.value = (-0.1).toFixed(4)
    if (z_input.value > 5) z_input.value = (5.0).toFixed(4)
}
z = parseFloat(z_input.value)

var line_restwave = line_data.data['restwave']
var waveshift_lines = (widgetinfos.data['waveframe_active'][0] == 0) ? 1+z : 1 ;
var waveshift_spec = (widgetinfos.data['waveframe_active'][0] == 0) ? 1 : 1/(1+z) ;

// Update x-range if in rest waveframe
if (widgetinfos.data['waveframe_active'][0] == 1) {
    // Trick: get previous z value from widgetinfos (which was not updated yet)
    var old_z = parseFloat(widgetinfos.data['z_input_value'][0]) ;
    var xrange_shift = (1+old_z)/(1+z) ;
    fig.x_range.start *= xrange_shift ;
    fig.x_range.end *= xrange_shift ;
}

// Update lines, overlap bands, data+noise wavelength array
shift_lines(lines, line_labels, line_restwave, waveshift_lines) ;
shift_lines(zlines, zline_labels, line_restwave, waveshift_lines) ;
if (overlap_bands.length>0) shift_bands(overlap_bands, waveshift_spec) ;
for(var i=0; i<spectra.length; i++) {
    shift_plotwave(spectra[i], waveshift_spec) ;
}
if (coaddcam_spec) shift_plotwave(coaddcam_spec, waveshift_spec) ;

// Update model wavelength array
// If othermodel is there, only the 'othermodel' curve is redshifted, and not the 'model' curve.
if (model) {
    var waveshift_model = (widgetinfos.data['waveframe_active'][0] == 0) ? 1 : 1/(1+z) ;
    if (othermodel === null) {
        var zmodel = 0.0 ;
        if(metadata.data['Z'] !== undefined) zmodel = metadata.data['Z'][ispectrumslider.value] ;
        waveshift_model *= ( (1+z)/(1+zmodel) ) ;
    }
    shift_plotwave(model, waveshift_model) ;
}
if (othermodel) {
    var zmodel = othermodel.data['zref'][0] ;
    var waveshift_model = (widgetinfos.data['waveframe_active'][0] == 0) ? (1+z)/(1+zmodel) : 1/(1+zmodel) ;
    shift_plotwave(othermodel, waveshift_model) ;
}

widgetinfos.data['z_input_value'][0] = z_input.value ;

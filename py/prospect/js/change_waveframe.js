// CustomJS, callback for the "Obs/rest waveframe" buttons
//  - Requires to include functions in shift_wave.js
//  - args = spectra, coaddcam_spec, model, othermodel, metadata, ispectrumslider,
//           z_input, widgetinfos, line_data, lines,
//           line_labels, zlines, zline_labels, overlap_waves, overlap_bands, fig.

var z = parseFloat(z_input.value)

// Change x-axis
if (cb_obj.active == 0) {
    fig.x_range.start = fig.x_range.start * (1+z) ;
    fig.x_range.end = fig.x_range.end * (1+z) ;
} else {
    fig.x_range.start = fig.x_range.start / (1+z) ;
    fig.x_range.end = fig.x_range.end / (1+z) ;
}

// Change x-position of glyphs
var line_restwave = line_data.data['restwave']
var waveshift_lines = (cb_obj.active == 0) ? 1+z : 1 ;
var waveshift_spec = (cb_obj.active == 0) ? 1 : 1/(1+z) ;

shift_lines(lines, line_labels, line_restwave, waveshift_lines) ;
shift_lines(zlines, zline_labels, line_restwave, waveshift_lines) ;

if (overlap_bands.length>0) shift_bands(overlap_bands, waveshift_spec) ;

for(var i=0; i<spectra.length; i++) {
    shift_plotwave(spectra[i], waveshift_spec) ;
}
if (coaddcam_spec) shift_plotwave(coaddcam_spec, waveshift_spec) ;

if (model) {
    var zmodel = 0.0 ;
    if(metadata.data['Z'] !== undefined) zmodel = metadata.data['Z'][ispectrumslider.value] ;
    var waveshift_model = (cb_obj.active == 0) ? (1+z)/(1+zmodel) : 1/(1+zmodel) ;
    shift_plotwave(model, waveshift_model) ;
}
if (othermodel) {
    var zmodel = othermodel.data['zref'][0] ;
    var waveshift_model = (cb_obj.active == 0) ? (1+z)/(1+zmodel) : 1/(1+zmodel) ;
    shift_plotwave(othermodel, waveshift_model) ;
}

widgetinfos.data['waveframe_active'][0] = cb_obj.active ;


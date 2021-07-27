// Utility functions used in CustomJS callbacks
//  Moving glyphs along the x-axis (change redshift or wavelength frame)


function shift_plotwave(cds_spec, waveshift) {
    var origwave = cds_spec.data['origwave']
    var plotwave = cds_spec.data['plotwave']
    if ( plotwave[0] != origwave[0] * waveshift ) { // Avoid redo calculation if not needed
        for (var j=0; j<plotwave.length; j++) {
            plotwave[j] = origwave[j] * waveshift ;
        }
        cds_spec.change.emit()
    }
}


function shift_lines(lines, line_labels, line_restwave, waveshift) {
    for (var i=0; i<line_restwave.length; i++) {
        lines[i].location = line_restwave[i] * waveshift ;
        line_labels[i].x = line_restwave[i] * waveshift ;
    }
}


function shift_bands(bands, band_waves, waveshift) {
    for (var i=0; i<bands.length; i++) {
        bands[i].left = band_waves[i][0] * waveshift ;
        bands[i].right = band_waves[i][1] * waveshift ;
    }
}


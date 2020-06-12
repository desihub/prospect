// Utility javascript function used in CustomJS callbacks
// Requires to include functions in: interp_grid.js

function coadd_brz_cameras(wave_in, flux_in, noise_in) {
    // Camera-coadd brz spectra.
    // each "_in" must have 3 entries (brz)
    // TODO handle case of no noise

    // Find b,r,z ordering in input arrays
    var wave_start = [wave_in[0][0], wave_in[1][0], wave_in[2][0]]
    var i_b = wave_start.indexOf(Math.min.apply(Math, wave_start))
    var i_z = wave_start.indexOf(Math.max.apply(Math, wave_start))
    var i_r = 1
    for (var i=0; i<3; i++) {
        if ( (i_b != i) && (i_z != i) ) i_r = i
    }

    var wave_out = []
    var flux_out = []
    var noise_out = []
    var margin = 20
    for (var i=0; i<wave_in[i_b].length; i++) { // b
        if (wave_in[i_b][i] < wave_in[i_b][wave_in[i_b].length-1] - margin) {
            wave_out.push(wave_in[i_b][i])
            flux_out.push(flux_in[i_b][i])
            noise_out.push(noise_in[i_b][i])
        }
    }
    var the_lim = wave_out[wave_out.length-1]
    for (var i=0; i<wave_in[i_r].length; i++) { // r
        if ( (wave_in[i_r][i] < wave_in[i_r][wave_in[i_r].length-1] - margin) && (wave_in[i_r][i] > the_lim)) {
            wave_out.push(wave_in[i_r][i])
            flux_out.push(flux_in[i_r][i])
            noise_out.push(noise_in[i_r][i])
        }
    }
    the_lim = wave_out[wave_out.length-1]
    for (var i=0; i<wave_in[i_z].length; i++) { // z
        if (wave_in[i_z][i] > the_lim) {
            wave_out.push(wave_in[i_z][i])
            flux_out.push(flux_in[i_z][i])
            noise_out.push(noise_in[i_z][i])
        }
    }
    for (var i=0; i<wave_out.length; i++) { // combine in overlapping regions
        var b1 = -1
        var b2 = -1
        if ( (wave_out[i] > wave_in[i_r][0]) && (wave_out[i] < wave_in[i_b][wave_in[i_b].length-1]) ) { // br
            b1 = 0
            b2 = 1
        }
        if ( (wave_out[i] > wave_in[i_z][0]) && (wave_out[i] < wave_in[i_r][wave_in[i_r].length-1]) ) {  // rz
            b1 = 1
            b2 = 2
        }
        if (b1 != -1) {
            var phi1 = interp_grid(wave_out[i], wave_in[b1], flux_in[b1])
            var noise1 = interp_grid(wave_out[i], wave_in[b1], noise_in[b1])
            var phi2 = interp_grid(wave_out[i], wave_in[b2], flux_in[b2])
            var noise2 = interp_grid(wave_out[i], wave_in[b2], noise_in[b2])
            if ( noise1 > 0 && noise2 > 0 ) {
                var iv1 = 1/(noise1*noise1)
                var iv2 = 1/(noise2*noise2)
                var iv = iv1+iv2
                noise_out[i] = 1/Math.sqrt(iv)
                flux_out[i] = (iv1*phi1+iv2*phi2)/iv
            }
        }
    }
    return [wave_out, flux_out, noise_out]
}

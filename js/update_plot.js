// update_plot() CustomJS


var d0 = new Date()
var t0 = d0.getTime()

var ifiber = ifiberslider.value
var nsmooth = smootherslider.value

// if (show_prev_vi_select.value == "Yes") {
//     vi_info_div.text = targetinfo.data['vi_info'][ifiber];
// } else {
//     vi_info_div.text = " ";
// }

if (cb_obj == ifiberslider) { // update VI widgets + infos for current spectrum

    // BYPASS DIV !!!
    // Code can clearly be better written... todo later...
//    target_info_div.text = targetinfo.data['target_info'][ifiber]
    targ_disp_cds.data['TARGETING'] = [ targetinfo.data['target_info'][ifiber] ]
    targ_disp_cds.data['mag_G'] = [ targetinfo.data['mag_G'][ifiber].toFixed(2) ]
    targ_disp_cds.data['mag_R'] = [ targetinfo.data['mag_R'][ifiber].toFixed(2) ]
    targ_disp_cds.data['mag_Z'] = [ targetinfo.data['mag_Z'][ifiber].toFixed(2) ]
    targ_disp_cds.data['mag_W1'] = [ targetinfo.data['mag_W1'][ifiber].toFixed(2) ]
    targ_disp_cds.data['mag_W2'] = [ targetinfo.data['mag_W2'][ifiber].toFixed(2) ]
    targ_disp_cds.change.emit()
    if(targetinfo.data['z'] != undefined && zcat_disp_cds != null) {
        zcat_disp_cds.data['SPECTYPE'] = [ targetinfo.data['spectype'][ifiber] ]
        zcat_disp_cds.data['Z'] = [ targetinfo.data['z'][ifiber].toFixed(4) ]
        zcat_disp_cds.data['ZERR'] = [ targetinfo.data['zerr'][ifiber].toFixed(4) ]
        zcat_disp_cds.data['ZWARN'] = [ targetinfo.data['zwarn'][ifiber] ]
        zcat_disp_cds.data['DeltaChi2'] = [ targetinfo.data['deltachi2'][ifiber].toFixed(1) ]
        zcat_disp_cds.change.emit()
    }
    
    vi_comment_input.value = targetinfo.data['VI_comment'][ifiber]
    vi_name_input.value = (targetinfo.data['VI_scanner'][ifiber]).trim()
    vi_class_input.active = vi_class_labels.indexOf(targetinfo.data['VI_class_flag'][ifiber]) // -1 if nothing
    var issues_on = []
    for (var i=0; i<vi_issue_slabels.length; i++) {
        if ( (targetinfo.data['VI_issue_flag'][ifiber]).indexOf(vi_issue_slabels[i]) >= 0 ) {
            issues_on.push(i)
        }
    }
    vi_issue_input.active = issues_on
    vi_z_input.value = (targetinfo.data['VI_z'][ifiber]).trim()
    vi_category_select.value = targetinfo.data['VI_spectype'][ifiber]
}

if(targetinfo.data['z'] != undefined && cb_obj == ifiberslider) {
    var z = targetinfo.data['z'][ifiber]
    var z1 = Math.floor(z*100) / 100
    zslider.value = z1
    dzslider.value = (z - z1)
}

function get_y_minmax(pmin, pmax, data) {
    // copy before sorting to not impact original, and filter out NaN
    var dx = data.slice().filter(Boolean)
    dx.sort()
    var imin = Math.floor(pmin * dx.length)
    var imax = Math.floor(pmax * dx.length)
    return [dx[imin], dx[imax]]
}

// Smoothing kernel
if (nsmooth > 0) {
    var kernel = [];
    for(var i=-2*nsmooth; i<=2*nsmooth; i++) {
        kernel.push(Math.exp(-(i**2)/(2*nsmooth)))
    }
    var kernel_offset = Math.floor(kernel.length/2)
}

function smooth_data(data_in, kernel, kernel_offset, quadrature=false) {
    // by default : out_j ~ (sum K_i in_i) / (sum K_i)
    // if quadrature is true (for noise) : out_j^2 ~ (sum K_i^2 in_i^2) / (sum K_i)^2
    var smoothed_data = data_in.slice()
    for (var j=0; j<data_in.length; j++) {
        smoothed_data[j] = 0.0
        var weight = 0.0
        // TODO: speed could be improved by moving `if` out of loop
        for (var k=0; k<kernel.length; k++) {
            var m = j+k-kernel_offset
            if((m >= 0) && (m < data_in.length)) {
                var fx = data_in[m]
                if(fx == fx) {
                    if (quadrature==true) {
                        smoothed_data[j] = smoothed_data[j] + (fx * kernel[k])**2
                    } else {
                        smoothed_data[j] = smoothed_data[j] + fx * kernel[k]
                    }
                    weight += kernel[k]
                }
            }
        }
        if (quadrature==true) {
            smoothed_data[j] = Math.sqrt(smoothed_data[j]) / weight
        } else {
            smoothed_data[j] = smoothed_data[j] / weight
        }
    }
    return smoothed_data
}

// Find nearest index in grid, left from point; use dichotomy method
function index_dichotomy(point, grid) {
    if ( point < grid[0] ) return 0
    if ( point > grid[grid.length-1] ) return grid.length-2
    var i_left = 0
    var i_center = 0
    var i_right = grid.length-1
    while ( i_right - i_left != 1) {
        i_center = i_left + Math.floor((i_right-i_left)/2)
        if ( point >= grid[i_center] ) {
            i_left = i_center
        } else {
            i_right = i_center
        }
    }
    return i_left
}

// Basic linear interpolation at on point
function interp_grid(xval, xarr, yarr) {
    var index = index_dichotomy(xval, xarr)
    var a = (yarr[index+1] - yarr[index])/(xarr[index+1] - xarr[index])
    var b = yarr[index]-a*xarr[index]
    var yval = a*xval+b
    return yval
}

// Coadd brz spectra. Similar to the python code mycoaddcam()
function coadd_brz_cams(wave_in, flux_in, noise_in) {
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

// Smooth plot and recalculate ymin/ymax
var ymin = 0.0
var ymax = 0.0
for (var i=0; i<spectra.length; i++) {
    var data = spectra[i].data
    var origflux = data['origflux'+ifiber]
    if ('plotnoise' in data) {
        var orignoise = data['orignoise'+ifiber]
    }
    if (nsmooth == 0) {
        data['plotflux'] = origflux.slice()
        if ('plotnoise' in data) {
            data['plotnoise'] = orignoise.slice()
        }
    } else {
        data['plotflux'] = smooth_data(origflux, kernel, kernel_offset)
        if ('plotnoise' in data) {
            // Add noise in quadrature
            data['plotnoise'] = smooth_data(orignoise, kernel, kernel_offset, quadrature=true)
        }
    }
    spectra[i].change.emit()

    tmp = get_y_minmax(0.01, 0.99, data['plotflux'])
    ymin = Math.min(ymin, tmp[0])
    ymax = Math.max(ymax, tmp[1])
}

// update camera-coadd
// Here I choose to do coaddition on the smoothed spectra (should be ok?)
if (coaddcam_spec) {
    var wave_in = []
    var flux_in = []
    var noise_in = []
    for (var i=0; i<3; i++) {
        var data = spectra[i].data
        wave_in.push(data['plotwave'].slice())
        flux_in.push(data['plotflux'].slice())
        if ('plotnoise' in data) {
            noise_in.push(data['plotnoise'].slice())
        } else {
            var dummy_noise = []
            for (var j=0; j<data['plotflux'].length; j++) dummy_noise.push(1)
            noise_in.push(dummy_noise)
        }
    }
    var coadd_infos = coadd_brz_cams(wave_in, flux_in, noise_in)
    coaddcam_spec.data['plotwave'] = coadd_infos[0].slice()
    coaddcam_spec.data['plotflux'] = coadd_infos[1].slice()
    coaddcam_spec.data['plotnoise'] = coadd_infos[2].slice()
    coaddcam_spec.change.emit()
}

// update model
if(model) {
    var origflux = model.data['origflux'+ifiber]
    if (nsmooth == 0) {
        model.data['plotflux'] = origflux.slice()
    } else {
        model.data['plotflux'] = smooth_data(origflux, kernel, kernel_offset)
    }
    model.change.emit()
}

// update y_range
if(ymin<0) {
    fig.y_range.start = ymin * 1.4
} else {
    fig.y_range.start = ymin * 0.6
}
fig.y_range.end = ymax * 1.4

// update target image
if (imfig_source) {
    imfig_source.data.url[0] = imfig_urls[ifiber][0]
    imfig_source.data.txt[0] = imfig_urls[ifiber][2]
    imfig_source.change.emit()
}


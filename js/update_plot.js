// update_plot() CustomJS
// Requires to include functions in: adapt_plotrange.js, coadd_brz_cameras.js, 
//     interp_grid.js, smooth_data.js

// TODO : optimize : if cd_obj==smootherslider don't need to do everything... (change imaging data etc)

var ifiber = ifiberslider.value
var nsmooth = smootherslider.value

if (cb_obj == ifiberslider) { // update VI widgets + infos for current spectrum

    // BYPASS DIV !!!
    // Code can clearly be better written... todo later...
//    target_info_div.text = targetinfo.data['target_info'][ifiber]
    targ_disp_cds.data['Target ID'] = [ targetinfo.data['targetid'][ifiber] ]
    targ_disp_cds.data['Target class'] = [ targetinfo.data['target_info'][ifiber] ]
    targ_disp_cds.data['mag_G'] = [ targetinfo.data['mag_G'][ifiber].toFixed(2) ]
    targ_disp_cds.data['mag_R'] = [ targetinfo.data['mag_R'][ifiber].toFixed(2) ]
    targ_disp_cds.data['mag_Z'] = [ targetinfo.data['mag_Z'][ifiber].toFixed(2) ]
    targ_disp_cds.data['mag_W1'] = [ targetinfo.data['mag_W1'][ifiber].toFixed(2) ]
    targ_disp_cds.data['mag_W2'] = [ targetinfo.data['mag_W2'][ifiber].toFixed(2) ]
    targ_disp_cds.change.emit()
    if(targetinfo.data['z'] != undefined && zcat_disp_cds != null) {
        if (fit_results != undefined) {
            zcat_disp_cds.data['SPECTYPE'] = fit_results['SPECTYPE'][ifiberslider.value].slice() // (0,num_best_fits)  
            zcat_disp_cds.data['SUBTYPE'] = fit_results['SUBTYPE'][ifiberslider.value].slice()
            zcat_disp_cds.data['Z'] = fit_results['Z'][ifiberslider.value].slice()
            zcat_disp_cds.data['ZERR'] = fit_results['ZERR'][ifiberslider.value].slice()
            zcat_disp_cds.data['ZWARN'] = fit_results['ZWARN'][ifiberslider.value].slice()
            var chi2s = fit_results['CHI2'][ifiberslider.value].slice() // Custom DeltaChi2 calculation
            var full_deltachi2s = []
            for (var i=0; i<fit_results['Nfit']-1; i++) {
                full_deltachi2s.push(chi2s[i+1]-chi2s[i])
            }
            full_deltachi2s.push(-1)
            zcat_disp_cds.data['DeltaChi2'] = full_deltachi2s
            for (var i=0; i<fit_results['Nfit']; i++) {
                zcat_disp_cds.data['Z'][i] = zcat_disp_cds.data['Z'][i].toFixed(4)
                zcat_disp_cds.data['ZERR'][i] = zcat_disp_cds.data['ZERR'][i].toFixed(4)
                zcat_disp_cds.data['DeltaChi2'][i] = zcat_disp_cds.data['DeltaChi2'][i].toFixed(1)
            }
        } else {
            zcat_disp_cds.data['SPECTYPE'] = [ targetinfo.data['spectype'][ifiber] ]
            zcat_disp_cds.data['SUBTYPE'] = [ targetinfo.data['subtype'][ifiber] ]
            zcat_disp_cds.data['Z'] = [ targetinfo.data['z'][ifiber].toFixed(4) ]
            zcat_disp_cds.data['ZERR'] = [ targetinfo.data['zerr'][ifiber].toFixed(4) ]
            zcat_disp_cds.data['ZWARN'] = [ targetinfo.data['zwarn'][ifiber] ]
            zcat_disp_cds.data['DeltaChi2'] = [ targetinfo.data['deltachi2'][ifiber].toFixed(1) ]
        }
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
    
    // update target image
    if (imfig_source) {
        imfig_source.data.url[0] = imfig_urls[ifiber][0]
        imfig_source.data.txt[0] = imfig_urls[ifiber][2]
        imfig_source.change.emit()
    }
}

if(targetinfo.data['z'] != undefined && cb_obj == ifiberslider && model_select == undefined) {
    // if model_select is defined : this will be done in select_model.
    var z = targetinfo.data['z'][ifiber]
    z_input.value = z.toFixed(4)
}

// Smoothing kernel
if (nsmooth > 0) {
    var kernel = [];
    for(var i=-2*nsmooth; i<=2*nsmooth; i++) {
        kernel.push(Math.exp(-(i**2)/(2*(nsmooth**2))))
    }
}

// Smooth plot and recalculate ymin/ymax TODO: ymin/max only if cb_obj==ifiberslider ?
var ymin = 0.0
var ymax = 0.0
var ivar_weight = true // switch for ivar-weighting in smoothing
for (var i=0; i<spectra.length; i++) {
    var data = spectra[i].data
    var origflux = data['origflux'+ifiber]
    var ivar = []
    if ('plotnoise' in data) {
        var orignoise = data['orignoise'+ifiber]
        for (var j=0; j<orignoise.length; j++) {
            if (orignoise[j]>0) {
                ivar.push( 1.0/(orignoise[j])**2 )
            } else {
                ivar.push(0)
            }
        }
    }
    if (nsmooth == 0) {
        data['plotflux'] = origflux.slice()
        if ('plotnoise' in data) {
            data['plotnoise'] = orignoise.slice()
        }
    } else {
        if ('plotnoise' in data) {
            data['plotflux'] = smooth_data(origflux, kernel, ivar_in=ivar, ivar_weight=ivar_weight)
            if (ivar_weight==true) {
                data['plotnoise'] = smooth_noise(ivar, kernel, ivar_weight=true)
            } else {
                data['plotnoise'] = smooth_noise(orignoise, kernel, ivar_weight=false)
            }
        } else {
            data['plotflux'] = smooth_data(origflux, kernel, ivar_weight=false)
        }
    }
    spectra[i].change.emit()

    var tmp = adapt_plotrange(0.01, 0.99, data['plotflux'])
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
    var coadd_infos = coadd_brz_cameras(wave_in, flux_in, noise_in)
    coaddcam_spec.data['plotwave'] = coadd_infos[0].slice()
    coaddcam_spec.data['plotflux'] = coadd_infos[1].slice()
    coaddcam_spec.data['plotnoise'] = coadd_infos[2].slice()
    coaddcam_spec.change.emit()
}

// update model
if (model) {
    var origflux = model.data['origflux'+ifiber]
    if (nsmooth == 0) {
        model.data['plotflux'] = origflux.slice()
    } else {
        model.data['plotflux'] = smooth_data(origflux, kernel)
    }
    model.change.emit()
}

// update other model
if (othermodel) {
    if (cb_obj == smootherslider) {
        var origflux = othermodel.data['origflux']
        if (nsmooth == 0) {
            othermodel.data['plotflux'] = origflux.slice()
        } else {
            othermodel.data['plotflux'] = smooth_data(origflux, kernel)
        }
        othermodel.change.emit()
    } else if (cb_obj == ifiberslider) {
        // Trick to trigger execution of select_model.js
        // Reset othermodel to best fit. Smoothing is done in select_model.js
        var trigger_value = model_select.options[0]
        if (model_select.value == trigger_value) {
            trigger_value = model_select.options[1]
        }
        model_select.value = trigger_value
        model_select.value = 'Best fit'
    }
}

// update y_range
if(ymin<0) {
    fig.y_range.start = ymin * 1.4
} else {
    fig.y_range.start = ymin * 0.6
}
fig.y_range.end = ymax * 1.4



// update_plot() CustomJS
// Requires to include functions in: adapt_plotrange.js, coadd_brz_cameras.js,
//     interp_grid.js, smooth_data.js
//
// TODO : optimize : if cd_obj==smootherslider don't need to do everything... (change imaging data etc)
//
// One of these will have changed.
//
var i_spectrum = ispectrumslider.value;
var nsmooth = smootherslider.value;

//
// Smoothing kernel.
//
var kernel = [];
if (nsmooth > 0)
    for(var i=-2*nsmooth; i<=2*nsmooth; i++)
        kernel.push(Math.exp(-(i**2)/(2*(nsmooth**2))));
//
// update VI widgets + infos for current spectrum
//
if (cb_obj == ispectrumslider) {
    //
    // Update metadata using "shortcds" objects.
    //
    var shortcds_list = [shortcds_table_a, shortcds_table_b, shortcds_table_c]
    if (shortcds_table_d !== null) {
        shortcds_list.push(shortcds_table_d)
    }
    //console.log(typeof(shortcds_table_d), shortcds_table_d != undefined, shortcds_table_d !== undefined, typeof(shortcds_table_d) !== undefined, shortcds_table_d != null, shortcds_table_d !== null)
    for (var i_cds=0; i_cds<shortcds_list.length; i_cds++) {
        var table_keys = Object.getOwnPropertyNames(shortcds_list[i_cds].data);
        for (var i=0; i<table_keys.length; i++) {
            if ( table_keys[i].includes('mag_') ) {
                shortcds_list[i_cds].data[table_keys[i]] = [ metadata.data[table_keys[i]][i_spectrum].toFixed(2) ];
            } else {
                shortcds_list[i_cds].data[table_keys[i]] = [ metadata.data[table_keys[i]][i_spectrum] ];
            }
        }
        shortcds_list[i_cds].change.emit();
    }
    if (shortcds_table_z !== null) {
        if (fit_results !== null) {
            shortcds_table_z.data['SPECTYPE'] = fit_results['SPECTYPE'][i_spectrum].slice(); // (0,num_best_fits)
            shortcds_table_z.data['SUBTYPE'] = fit_results['SUBTYPE'][i_spectrum].slice();
            shortcds_table_z.data['Z'] = fit_results['Z'][i_spectrum].slice();
            shortcds_table_z.data['ZERR'] = fit_results['ZERR'][i_spectrum].slice();
            shortcds_table_z.data['ZWARN'] = fit_results['ZWARN'][i_spectrum].slice();
            var chi2s = fit_results['CHI2'][i_spectrum].slice(); // Custom DeltaChi2 calculation
            var full_deltachi2s = [];
            for (var i=0; i<fit_results['Nfit']-1; i++) {
                full_deltachi2s.push(chi2s[i+1]-chi2s[i]);
            }
            full_deltachi2s.push(-1);
            shortcds_table_z.data['DELTACHI2'] = full_deltachi2s;
            for (var i=0; i<fit_results['Nfit']; i++) {
                shortcds_table_z.data['Z'][i] = shortcds_table_z.data['Z'][i].toFixed(4);
                shortcds_table_z.data['ZERR'][i] = shortcds_table_z.data['ZERR'][i].toFixed(4);
                shortcds_table_z.data['DELTACHI2'][i] = shortcds_table_z.data['DELTACHI2'][i].toFixed(1);
            }
        } else {
            var table_keys = Object.getOwnPropertyNames(shortcds_table_z.data);
            for (var i=0; i<table_keys.length; i++) {
                if ( table_keys[i].includes('CHI2') ) {
                    shortcds_table_z.data[table_keys[i]] = [ metadata.data[table_keys[i]][i_spectrum].toFixed(1) ];
                } else if ( table_keys[i].includes('CLASS') || table_keys[i].includes('TYPE') || table_keys[i].includes('WARN')) {
                    shortcds_table_z.data[table_keys[i]] = [ metadata.data[table_keys[i]][i_spectrum] ];
                } else { // Z, ZERR
                    shortcds_table_z.data[table_keys[i]] = [ metadata.data[table_keys[i]][i_spectrum].toFixed(4) ];
                }
            }
        }
        shortcds_table_z.change.emit();
    }
    //
    // Update VI
    //
    vi_std_comment_select.value = ' ';
    vi_comment_input.value = metadata.data['VI_comment'][i_spectrum];
    vi_name_input.value = (metadata.data['VI_scanner'][i_spectrum]).trim();
    vi_quality_input.active = vi_quality_labels.indexOf(metadata.data['VI_quality_flag'][i_spectrum]); // -1 if nothing
    var issues_on = [];
    for (var i=0; i<vi_issue_slabels.length; i++) {
        if ( (metadata.data['VI_issue_flag'][i_spectrum]).indexOf(vi_issue_slabels[i]) >= 0 ) {
            issues_on.push(i);
        }
    }
    vi_issue_input.active = issues_on;
    vi_z_input.value = (metadata.data['VI_z'][i_spectrum]).trim();
    if (metadata.data['VI_spectype'][i_spectrum] == '' ) {
        vi_category_select.value = ' '; // setting value='' does not work well with Select
    } else {
        vi_category_select.value = metadata.data['VI_spectype'][i_spectrum];
    }

    //
    // update target image
    //
    if (imfig_source) {
        imfig_source.data.url[0] = imfig_urls[i_spectrum][0];
        imfig_source.data.txt[0] = imfig_urls[i_spectrum][2];
        imfig_source.change.emit();
    }
    //
    // reset x-range
    //
    fig.x_range.start = xrange[0];
    fig.x_range.end = xrange[1];
    if (widgetinfos.data['waveframe_active'][0] == 1) {
        // Combines with the x-range update done in change_redshift.js
        var old_z = parseFloat(z_input.value)
        fig.x_range.start /= (1+old_z);
        fig.x_range.end /= (1+old_z);
    }
}
//
// Update redshift.
//
if (metadata.data['Z'] !== undefined && cb_obj == ispectrumslider && model_select === null) {
    // if model_select is defined : this will be done in select_model.
    z_input.value = metadata.data['Z'][i_spectrum].toFixed(4);
}
//
// Smooth plot and recalculate ymin/ymax.
//
var ymin = 0.0;
var ymax = 0.0;
// Switch for ivar-weighting in smoothing. We normally want this true unless
// nsmooth == 0, to simplify some of the logic below.
var ivar_weight = nsmooth > 0;
var valid_y_range = false;
for (var i=0; i<spectra.length; i++) {
    // console.log("Updating spectrum " + (i+1) + " (of " + spectra.length + ") of object " + i_spectrum + ".");
    var data = spectra[i].data;
    var origflux = data['origflux'+String(i_spectrum)];
    if (origflux.filter(isFinite).length == 0) {
        alert("Spectrum " + (i+1) + " (of " + spectra.length + ") of object " + i_spectrum + " has no valid data!");
        data["plotflux"] = (function(){ var foo = []; for (var j=0; j<origflux.length; j++) foo.push(1.0); return foo;})();
        if ("plotnoise" in data) data["plotnoise"] = data["plotflux"].slice();
    } else {
        if ('plotnoise' in data) {
            var orignoise = data['orignoise'+String(i_spectrum)];
            var ivar = [];
            for (var j=0; j<orignoise.length; j++) ivar.push( orignoise[j] > 0 ? 1.0/(orignoise[j])**2 : 0 );
            data["plotflux"] = smooth_data(origflux, kernel, {"ivar_in": ivar, "ivar_weight": ivar_weight});
            data["plotnoise"] = smooth_noise((ivar_weight ? ivar : orignoise), kernel, {"ivar_weight": ivar_weight});
        } else {
            data["plotflux"] = smooth_data(origflux, kernel, {"ivar_weight": false});
        }
        var y_vect = []; // ymin/ymax calculation is based on pixels in fig.x_range
        for (var i_pix=0; i_pix<data['plotflux'].length; i_pix++) {
            if ( (data['plotwave'][i_pix] > fig.x_range.start ) && (data['plotwave'][i_pix] < fig.x_range.end) ) {
                y_vect.push(data['plotflux'][i_pix])
            }
        }
        if (y_vect.length > 0) {
            var tmp = adapt_plotrange(0.01, 0.99, y_vect);
            ymin = Math.min(ymin, tmp[0]);
            ymax = Math.max(ymax, tmp[1]);
        }
        valid_y_range = isFinite(ymin) && isFinite(ymax);
    }
    spectra[i].change.emit();
}
//
// After initial limits are calculated by adapt_plotrange(), add some
// extra padding.
//
if (valid_y_range) {
    fig.y_range.start = ymin < 0 ? ymin * 1.4 : ymin * 0.6;
    fig.y_range.end = ymax * 1.4;
}
//
// update camera-coadd
// Here I choose to do coaddition on the smoothed spectra (should be ok?)
//
if (coaddcam_spec) {
    var wave_in = [];
    var flux_in = [];
    var noise_in = [];
    for (var i=0; i<3; i++) {
        var data = spectra[i].data;
        wave_in.push(data['plotwave'].slice());
        flux_in.push(data['plotflux'].slice());
        if ('plotnoise' in data) {
            noise_in.push(data['plotnoise'].slice());
        } else {
            var dummy_noise = [];
            for (var j=0; j<data['plotflux'].length; j++) dummy_noise.push(1);
            noise_in.push(dummy_noise);
        }
    }
    var coadd_infos = coadd_brz_cameras(wave_in, flux_in, noise_in);
    coaddcam_spec.data['plotwave'] = coadd_infos[0].slice();
    coaddcam_spec.data['plotflux'] = coadd_infos[1].slice();
    coaddcam_spec.data['plotnoise'] = coadd_infos[2].slice();
    coaddcam_spec.change.emit();
}
//
// update model
//
if (model) {
    model.data['plotflux'] = smooth_data(model.data['origflux'+String(i_spectrum)], kernel, {});
    model.change.emit();
}
//
// update other model
//
if (othermodel) {
    if (cb_obj == smootherslider) {
        othermodel.data['plotflux'] = smooth_data(othermodel.data['origflux'], kernel, {});
        othermodel.change.emit();
    } else if (cb_obj == ispectrumslider) {
        // Trick to trigger execution of select_model.js
        // Reset othermodel to best fit. Smoothing is done in select_model.js
        var trigger_value = model_select.options[0];
        if (model_select.value == trigger_value) {
            trigger_value = model_select.options[1];
        }
        model_select.value = trigger_value;
        model_select.value = 'Best fit';
    }
}

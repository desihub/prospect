// CustomJS, callback for the model selection
//  - Requires to include function in: interp_grid.js, smooth_data.js
//  - args = ispectrumslider, model_select, fit_templates, cds_model_2ndfit, cds_model, z_input
//             cds_othermodel, rrdetails, std_templates, median_spectra, smootherslider,
//             cds_metadata
// IN DEVLPT
// values for model_select (model_options) are hardcoded
// TODO : add smootherslider info (nsmooth, see code in update_plot )


function median(numbers) {
    var median = 0
    var arr = numbers.slice()
    arr.sort()
    if ( arr.length % 2 === 0 ) { // is even
        median = (arr[arr.length / 2 - 1] + arr[arr.length / 2]) / 2
    } else { // is odd
        median = arr[(arr.length - 1) / 2]
    }
    return median
}

var i_spectrum = ispectrumslider.value
var spec_z = 0.0

if (model_select.value == 'Best fit') {
    var origflux = cds_model.data['origflux'+String(i_spectrum)]

    cds_othermodel.data['origflux'] = origflux.slice()
    cds_othermodel.data['plotflux'] = origflux.slice()
    cds_othermodel.data['plotwave'] = cds_model.data['plotwave'].slice()
    cds_othermodel.data['origwave'] = cds_model.data['origwave'].slice()
    spec_z = cds_metadata.data['Z'][i_spectrum]

} else if (model_select.value == '2nd best fit') {
    var origflux = cds_model_2ndfit.data['origflux'+String(i_spectrum)]
    cds_othermodel.data['origflux'] = origflux.slice()
    cds_othermodel.data['plotflux'] = origflux.slice()
    cds_othermodel.data['plotwave'] = cds_model_2ndfit.data['plotwave']
    cds_othermodel.data['origwave'] = cds_model_2ndfit.data['origwave']
    spec_z = rrdetails['Z'][i_spectrum][1] // entry "1" => 2nd best fit

} else if ( (model_select.value).search("STD ") == 0) {
    // TODO : adapt median to relevant waverange (fct of z_spec at least).

    // Option (A): plot std template at an arbitrary/fixed redshift
    //if ( (model_select.value).search("GALAXY") > -1) spec_z = 0.7
    //if ( (model_select.value).search("QSO") > -1) spec_z = 1.5
    // Option (B): plot std template at the redshift given by z_input
    spec_z = parseFloat(z_input.value)

    var template_key = (model_select.value).slice(4)
    var shifted_template_wave = std_templates["wave_"+template_key].slice()
    // Median is computed in waverange [3600-9800] angstrom (RoI)
    var lambda_min_roi = 3600
    var lambda_max_roi = 9800
    var template_in_roi = []
    for (var i=0; i<shifted_template_wave.length; i++) {
        shifted_template_wave[i] = shifted_template_wave[i]*(1+spec_z)
        if (shifted_template_wave[i]>lambda_min_roi && shifted_template_wave[i]<lambda_max_roi) {
            template_in_roi.push(std_templates["flux_"+template_key][i])
        }
    }

    var median_goal = median_spectra.data['median'][i_spectrum]
    var median_template = median(template_in_roi)
    var scaled_template_flux = std_templates["flux_"+template_key].slice()
    for (var i=0; i<scaled_template_flux.length; i++) {
        scaled_template_flux[i] = scaled_template_flux[i]*median_goal/median_template
    }
    cds_othermodel.data['plotflux'] = scaled_template_flux
    cds_othermodel.data['origflux'] = scaled_template_flux.slice()

    cds_othermodel.data['plotwave'] = shifted_template_wave.slice()
    cds_othermodel.data['origwave'] = shifted_template_wave.slice()

} else { // recompute Nth best fit
    var i_fit = parseInt(model_select.value[0])-1 // hardcoded ("1st fit" => "1" => entry 0. Assumes N<10)
    var coeffs = rrdetails['COEFF'][i_spectrum][i_fit]
    var spectype = rrdetails['SPECTYPE'][i_spectrum][i_fit]
    var subtype = rrdetails['SUBTYPE'][i_spectrum][i_fit]
    spec_z = rrdetails['Z'][i_spectrum][i_fit]
    var template_wave = fit_templates["wave_"+spectype+"_"+subtype]
    var template_flux = fit_templates["flux_"+spectype+"_"+subtype]
    
    var model_wave = template_wave.slice()
    var model_flux = template_wave.slice()
    for (var j=0; j<model_wave.length; j++) {
        model_wave[j] = model_wave[j]*(1+spec_z)
        model_flux[j] = 0
        for (var i=0; i<template_flux.length; i++) {
            model_flux[j] += (coeffs[i]*template_flux[i][j])
        }
    }

    cds_othermodel.data['plotwave'] = model_wave.slice()
    cds_othermodel.data['origwave'] = model_wave.slice()
    cds_othermodel.data['origflux'] = model_flux.slice()
    cds_othermodel.data['plotflux'] = model_flux.slice()
}

var zref_vect = (cds_othermodel.data['plotwave']).slice() // trick to keep track of spec_z (see change_redshift.js)
for (var j=0; j<zref_vect.length; j++) zref_vect[j] = spec_z
cds_othermodel.data['zref'] = zref_vect.slice()

// Smooth plotflux
// Take into account binning difference between data and othermodel:
var binning_data = 0.8 ;
if (cds_model !== null) {
    binning_data = cds_model.data['plotwave'][1]-cds_model.data['plotwave'][0] ;
}
var binning_othermodel = cds_othermodel.data['plotwave'][1]-cds_othermodel.data['plotwave'][0] ;
var nsmooth = smootherslider.value * binning_data / binning_othermodel ;
if (nsmooth > 0) {
    var kernel = get_kernel(nsmooth);
    cds_othermodel.data['plotflux'] = smooth_data(cds_othermodel.data['origflux'], kernel, {})
}

// Change value for z_input
z_input.value = spec_z.toFixed(4)

cds_othermodel.change.emit()

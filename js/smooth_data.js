// Utility javascript function used in CustomJS callbacks
//  (Standalone module)

function smooth_data(data_in, kernel, ivar_in=0, ivar_weight=false, kernel_offset="default") {
    // by default: out_j ~ (sum K_i in_i) / (sum K_i)
    // ivar-weighting: out_j ~ (sum K_i in_i ivar_i) / (sum K_i ivar_i) ~ smooth(in*ivar)/smooth(ivar)
    // default kernel offset: symetric kernel around 0
    if (kernel_offset=="default") {
        kernel_offset = Math.floor(kernel.length/2)
    }
    var smoothed_data = data_in.slice()
    for (var j=0; j<data_in.length; j++) {
        smoothed_data[j] = 0.0
        var weight = 0.0
        // TODO: speed could be improved by moving `if` out of loop
        for (var k=0; k<kernel.length; k++) {
            var m = j+k-kernel_offset
            if((m >= 0) && (m < data_in.length)) {
                if (ivar_weight==false) {
                    var fx = data_in[m]
                    if(fx == fx) {
                        smoothed_data[j] = smoothed_data[j] + fx * kernel[k]
                        weight += kernel[k]
                    }
                } else {
                    var fx = data_in[m] * ivar_in[m]
                    if(fx == fx) {
                        smoothed_data[j] = smoothed_data[j] + fx * kernel[k]
                        weight += kernel[k]*ivar_in[m]
                    }
                }
            }
        }
        smoothed_data[j] = smoothed_data[j] / weight
    }
    return smoothed_data
}

function smooth_noise(noise_in, kernel, ivar_weight=false, kernel_offset="default") {
    // Add noise in quadrature : 
    // by default: out_j^2 ~ (sum K_i^2 in_i^2) / (sum K_i)^2
    // ivar-weighting: then noise_in is understood as ivar, and out_j ~ (sum K_i^2 in_i) / (sum K_i in_i)^2
    if (kernel_offset=="default") {
        kernel_offset = Math.floor(kernel.length/2)
    }
    var smoothed_noise = noise_in.slice()
    for (var j=0; j<noise_in.length; j++) {
        smoothed_noise[j] = 0.0
        var weight = 0.0
        for (var k=0; k<kernel.length; k++) {
            var m = j+k-kernel_offset
            if((m >= 0) && (m < noise_in.length)) {
                if (ivar_weight==false) {
                    var fx = noise_in[m]
                    if(fx == fx) {
                        smoothed_noise[j] = smoothed_noise[j] + (fx * kernel[k])**2
                        weight += kernel[k]
                    }
                } else {
                    var fx = kernel[k] * noise_in[m]
                    if(fx == fx) {
                        smoothed_noise[j] = smoothed_noise[j] + (fx * kernel[k])
                        weight += fx
                    }
                }
            }
        }
        smoothed_noise[j] = Math.sqrt(smoothed_noise[j]) / weight
    }
    return smoothed_noise
}

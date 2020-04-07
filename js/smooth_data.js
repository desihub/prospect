// Utility javascript function used in CustomJS callbacks
//  (Standalone module)

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

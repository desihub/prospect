// Utility javascript functions used in CustomJS callbacks
//  (Standalone module)

function index_dichotomy(point, grid) {
    // Find nearest index in grid, left from point; use dichotomy method
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

function interp_grid(xval, xarr, yarr) {
    // Basic linear interpolation of [xarr,yarr] on point xval
    var index = index_dichotomy(xval, xarr)
    var a = (yarr[index+1] - yarr[index])/(xarr[index+1] - xarr[index])
    var b = yarr[index]-a*xarr[index]
    var yval = a*xval+b
    return yval
}

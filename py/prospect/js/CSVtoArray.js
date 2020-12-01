// Utility javascript function
//  (Standalone module)

// Return array of string values
function CSVtoArray(text) {        
    var re_value = /(?!\s*$)\s*(?:'([^'\\]*(?:\\[\S\s][^'\\]*)*)'|"([^"\\]*(?:\\[\S\s][^"\\]*)*)"|([^,'"\s\\]*(?:\s+[^,'"\s\\]+)*))\s*(?:,|$)/g
    var a = []
    text.replace(re_value, // "Walk" the string using replace with callback.
        function(m0, m1, m2, m3) {
            // Remove backslash from \' in single quoted values.
            if      (m1 !== undefined) a.push(m1.replace(/\\'/g, "'"))
            // Remove backslash from \" in double quoted values.
            else if (m2 !== undefined) a.push(m2.replace(/\\"/g, '"'))
            else if (m3 !== undefined) a.push(m3)
            return '' // Return empty string.
        })
    // Handle special case of empty last value.
    if (/,\s*$/.test(text)) a.push('')
    return a
}

// CustomJS, storing VI information to browser's localstorage


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

function autosave_vi(title, vi_file_fields, cds_data) {
    // title : localStorage key where VI info is stored
    // vi_file_fields : VI fields to be stored (as defined in utils_specviewer)
    // cds_data : data from Bokeh CDS containing VI informations
    //            must contain at least "VI_class_flag", "VI_comment", "VI_issue_flag"
    var nb_fields = vi_file_fields.length
    var nspec = cds_data['VI_class_flag'].length
    
    var array_to_store = []
    // Get previous localStorage beforehand (avoid overwritting)
    var previously_stored_ispecs = []
    if (title in localStorage) {
        var recovered_csv = localStorage.getItem(title)
        var recovered_entries = recovered_csv.split("\n")
        for (var j=0; j<recovered_entries.length; j++) {
            var row = CSVtoArray(recovered_entries[j])
            array_to_store.push(row)
            previously_stored_ispecs.push(Number(row[0]))
        }
    }
    for (var i_spec=0; i_spec<nspec; i_spec++) {
         // Choice : Record only information if a VI classification was assigned
        if ( (cds_data['VI_class_flag'][i_spec] != "-1") ||
            (cds_data['VI_comment'][i_spec].trim() != "") ||
            (cds_data['VI_issue_flag'][i_spec].trim() != "") ||
            (cds_data['VI_z'][i_spec].trim() != "") ) {
            var row = [ i_spec.toString() ] // Store also spectrum number
            for (var j=0; j<vi_file_fields.length; j++) {
                var entry = cds_data[vi_file_fields[j][1]][i_spec]
                if ( typeof(entry)!="string" ) entry = entry.toString()
                entry = entry.replace(/"/g, '""')
                entry = entry.replace(/,/g, '","')
                row.push(entry)
            }
            i_rec = previously_stored_ispecs.indexOf(i_spec)
            if (i_rec == -1) {
                array_to_store.push(row)
            } else {
                array_to_store[i_rec] = row // Replace entry
            }
        }
    }

    var csv_to_store = ''
    for (var j=0; j<array_to_store.length; j++) {
        var row = (array_to_store[j]).join(' , ')
        csv_to_store += ( row.concat("\n") )
    }

    if (typeof(localStorage) !== "undefined") {
        localStorage.setItem(title, csv_to_store)
    } else {
        console.log("Warning : no local storage available in browser.")
    }
}


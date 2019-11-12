// CustomJS, storing VI information to browser's localstorage

function autosave_vi(title, vi_file_fields, cds_data) {
    // title : localStorage key where VI info is stored
    // vi_file_fields : VI fields to be stored (as defined in utils_specviewer)
    // cds_data : data from Bokeh CDS containing VI informations
    //            must contain at least "VI_class_flag", "VI_comment", "VI_issue_flag"
    var nb_fields = vi_file_fields.length
    var nspec = cds_data['VI_class_flag'].length
    
    var array_to_store = []        
    for (var i_spec=0; i_spec<nspec; i_spec++) {
         // Record only information if a VI classification was assigned
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
            array_to_store.push(row)
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


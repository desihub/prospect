// Functions used in several CustomJS callbacks
//  storing VI information to either localStorage or csv file
// Requires to include functions in: CSVtoArray.js

// Create a string variable containing VI records in CSV format.
function vi_to_csv(vi_file_fields, cds_data, for_localStorage, localStorage_key) {
    // vi_file_fields: VI fields to be stored (defined in prospect.utilities)
    // cds_data: data from Bokeh CDS containing VI informations
    //    must contain at least "VI_class_flag", "VI_comment", "VI_issue_flag"
    // for_localStorage (bool): if true, the output format is slightly modified:
    //    no header, add a first column providing spectrum number.
    // localStorage_key: if not undefined, previous VI information from the
    //    browser's localStorage (with corresponding key) is read
    //    beforehand and included to the output.

    var nb_fields = vi_file_fields.length
    var nspec = cds_data['VI_class_flag'].length

    var array_to_store = []

    if (for_localStorage == false) {
        var header = []
        for (var j=0; j<nb_fields; j++) header.push(vi_file_fields[j][0])
        array_to_store.push(header)
    }

    if (localStorage_key != undefined) {
        var previously_stored_ispecs = []
        if (localStorage_key in localStorage) {
            var recovered_csv = localStorage.getItem(localStorage_key)
            var recovered_entries = recovered_csv.split("\n")
            for (var j=0; j<recovered_entries.length; j++) {
                var row = CSVtoArray(recovered_entries[j])
                array_to_store.push(row)
                previously_stored_ispecs.push(Number(row[0]))
            }
        }
    }

    for (var i_spec=0; i_spec<nspec; i_spec++) {
        //  Record only information if a VI classification was assigned
        //    or some VI comment/issue/z was given:
        if ( (cds_data['VI_class_flag'][i_spec] != "-1") ||
            (cds_data['VI_comment'][i_spec].trim() != "") ||
            (cds_data['VI_issue_flag'][i_spec].trim() != "") ||
            (cds_data['VI_z'][i_spec].trim() != "") ) {
            var row = []
            if (for_localStorage == true) {
                row.push(i_spec.toString())
            }
            for (var j=0; j<vi_file_fields.length; j++) {
                var entry = cds_data[vi_file_fields[j][1]][i_spec]
                if (vi_file_fields[j][1] == "z") entry = entry.toFixed(4)
                if ( typeof(entry) != "string" ) entry = entry.toString()
                entry = entry.replace(/"/g, '""')
                entry = entry.replace(/,/g, '","')
                if (for_localStorage == false) {
                    if (entry==" ") entry = ""
                }
                row.push(entry)
            }
            var i_rec = -1
            if (localStorage_key != undefined) {
                var i_rec = -1
                i_rec = previously_stored_ispecs.indexOf(i_spec)
                if (i_rec == -1) {
                    array_to_store.push(row)
                } else {
                    array_to_store[i_rec] = row // Replace entry
                }
            } else {
                array_to_store.push(row)
            }
        }
    }

    var csv_to_store = ''
    for (var j=0; j<array_to_store.length; j++) {
        var row = (array_to_store[j]).join(',')
        csv_to_store += ( row.concat("\n") )
    }

    return csv_to_store
}

// Store VI information to localStorage
function autosave_vi_localStorage(vi_file_fields, cds_data, localStorage_key) {
    if (typeof(localStorage) !== "undefined") {
        var for_localStorage = true
        var csv_to_store = vi_to_csv(vi_file_fields, cds_data, for_localStorage, localStorage_key)
        localStorage.setItem(localStorage_key, csv_to_store)
    } else {
        console.log("Warning : no local storage available in browser.")
    }
}

// Store VI information to VI csv file
function download_vi_file(vi_file_fields, cds_data, output_file) {
    var for_localStorage = false
    var csv_to_store = vi_to_csv(vi_file_fields, cds_data, for_localStorage)
    var blob = new window.Blob([csv_to_store], {type: 'text/csv'})
    saveAs(blob, output_file)
}

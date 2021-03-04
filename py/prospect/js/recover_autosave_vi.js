// recover_vi_callback() CustomJS
// Requires to include functions in: CSVtoArray.js

// Recover auto-saved VI infos from browser's localStorage
// args = title, cds_metadata, vi_file_fields, ifiber, 
//   vi_comment_input, vi_name_input, vi_class_input, vi_issue_input, vi_issue_slabels, vi_class_labels

if (title in localStorage) {
    var recovered_csv = localStorage.getItem(title)
    var recovered_entries = recovered_csv.split("\n")
    for (var j=0; j<recovered_entries.length; j++) {
        var row = CSVtoArray(recovered_entries[j])
        var i_spec = Number(row[0])
        for (var k=1; k<row.length; k++) {
            if (vi_file_fields[k-1][1].includes('VI')) {
                cds_metadata.data[vi_file_fields[k-1][1]][i_spec] = row[k]
            }
        }
        if (i_spec == ifiber) { // update VI buttons for current spectrum
            vi_comment_input.value = cds_metadata.data['VI_comment'][ifiber] ;
            vi_name_input.value = cds_metadata.data['VI_scanner'][ifiber] ;
            vi_class_input.active = vi_class_labels.indexOf(cds_metadata.data['VI_class_flag'][ifiber]) ; // -1 if nothing
            var issues_on = []
            for (var i=0; i<vi_issue_slabels.length; i++) {
                if ( (cds_metadata.data['VI_issue_flag'][ifiber]).indexOf(vi_issue_slabels[i]) >= 0 ) {
                    issues_on.push(i)
                }
            }
            vi_issue_input.active = issues_on
        }
    }
    cds_metadata.change.emit()
}

// recover_vi_callback() CustomJS
// Requires to include functions in: CSVtoArray.js

// Recover auto-saved VI infos from browser's localStorage
// args = title, cds_metadata, output_file_fields, ispectrumslider,
//   vi_comment_input, vi_name_input, vi_quality_input, vi_issue_input, vi_issue_slabels, vi_quality_labels

if (title in localStorage) {
    var recovered_csv = localStorage.getItem(title)
    var recovered_entries = recovered_csv.split("\n")
    for (var j=0; j<recovered_entries.length; j++) {
        var row = CSVtoArray(recovered_entries[j])
        var i_spec = Number(row[0])
        for (var k=1; k<row.length; k++) {
            if (output_file_fields[k-1][1].includes('VI')) {
                cds_metadata.data[output_file_fields[k-1][1]][i_spec] = row[k]
            }
        }
        if (i_spec == ispectrumslider.value) { // update VI buttons for current spectrum
            vi_comment_input.value = cds_metadata.data['VI_comment'][i_spec] ;
            vi_name_input.value = cds_metadata.data['VI_scanner'][i_spec] ;
            vi_quality_input.active = vi_quality_labels.indexOf(cds_metadata.data['VI_quality_flag'][i_spec]) ; // -1 if nothing
            var issues_on = []
            for (var i=0; i<vi_issue_slabels.length; i++) {
                if ( (cds_metadata.data['VI_issue_flag'][i_spec]).indexOf(vi_issue_slabels[i]) >= 0 ) {
                    issues_on.push(i)
                }
            }
            vi_issue_input.active = issues_on
        }
    }
    cds_metadata.change.emit()
}

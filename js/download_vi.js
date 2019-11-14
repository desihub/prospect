// CustomJS, download VI results to file
//  args = cds_targetinfo, vi_file_fields, vi_filename

var nb_fields = vi_file_fields.length
var nspec = cds_targetinfo.data['VI_class_flag'].length

var array_to_store = []        
var header = []
for (var j=0; j<nb_fields; j++) header.push(vi_file_fields[j][0])
array_to_store.push(header)
for (var i=0; i<nspec; i++) {
     // Record only information if a VI classification was assigned
    if ( (cds_targetinfo.data['VI_class_flag'][i] != "-1") ||
         (cds_targetinfo.data['VI_comment'][i].trim() != "") ||
         (cds_targetinfo.data['VI_issue_flag'][i].trim() != "" ) ) {
        var row = []
        for (var j=0; j<vi_file_fields.length; j++) {
            var entry = cds_targetinfo.data[vi_file_fields[j][1]][i]
            if (vi_file_fields[j][1]=="z") entry = entry.toFixed(3)
            if ( typeof(entry)!="string" ) entry = entry.toString()
            entry = entry.replace(/"/g, '""')
            entry = entry.replace(/,/g, '","')
            if (entry=="" || entry==" ") entry = "--"
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

var blob = new window.Blob([csv_to_store], {type: 'text/csv'})
saveAs(blob, vi_filename)

// Old function replaced by saveAs : download(filename, csv_to_store)
//         function download(filename, text) {
//             var element = document.createElement('a')
//             element.setAttribute('href', 'data:text/csv;charset=utf-8,' + encodeURIComponent(text))
//             element.setAttribute('download', filename)
//             element.style.display = 'none'
//             document.body.appendChild(element)
//             element.click()
//             document.body.removeChild(element)
//         }

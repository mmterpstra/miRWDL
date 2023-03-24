version 1.0

task MultiQC {
    input {
        #should be files or dirs
        Array[File] files
        #Scales badly
        Int? memoryGb = "4"
	    String? multiQCModule  = "MultiQC/1.9-foss-2020a-Python-3.8.2"
        Int timeMinutes = 50
    }

    command {
        set -e
        module load ${multiQCModule} && multiqc --force --file-list ${write_lines(files)}
    }

    output {
        File dir =  "multiqc_data"
        File html = "multiqc_report.html"
    }

    runtime {
        memory: select_first([memoryGb * 1024,1024])
        timeMinutes: timeMinutes
    }
}
version 1.0

task FastQC {
    input {
        File inputFastq
        Int? memoryGb = "1"
	    String? fastqcModule  = "FastQC/0.11.9-Java-11"
        Int timeMinutes = 1 + ceil(size(inputFastq, "G")) * 10
    }

    command {
        set -e
        module load ${fastqcModule} && fastqc -o ./ ${inputFastq}
    }

    output {
        String fastqcDir =  basename(inputFastq, ".fastq.gz")
        #ile out_html1 = basename(inputFastq1, ".fastq.gz") + "_fastqc.html"
        File outZip = basename(inputFastq, ".fastq.gz") + "_fastqc.zip"
    }

    runtime {
        memory: select_first([memoryGb * 1024,1024])
        timeMinutes: timeMinutes
    }
}
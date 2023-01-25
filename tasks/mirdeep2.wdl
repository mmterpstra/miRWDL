version 1.0

task Quantifier {
    input {
        File inputCollapsedFasta
        Int? memoryGb = 1 + ceil(size(inputCollapsedFasta, "G"))*1
	    String mirdeepModule  = "mirdeep2/0.1.3-GCC-10.2.0-Perl-5.32.0"
        Int timeMinutes = 1 + ceil(size(inputFastq, "G")) * 10
    }

    command {
        set -e
        module load ${mirdeepModule} && quantifier.pl ${inputCollapsedFasta}
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

task QuantifierHumanDb {
    input {
        File inputCollapsedFasta
        Int? memoryGb = 1 + ceil(size(inputCollapsedFasta, "G"))*1
	    String miRDeepModule  = "mirdeep2/0.1.3-GCC-10.2.0-Perl-5.32.0"
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
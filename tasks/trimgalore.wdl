version 1.0

task TrimGalore {
    input {
        File inputFastq1
        File? inputFastq2
        File outputFastq1
        File? outputFastq2
        Int? memoryGb = 1

        #File? fastq_input_umi
        #String samplename
        #String identifier
	    String trimgaloreModule = "Trim_Galore/0.6.6-GCCcore-9.3.0-Python-3.8.2"
        #String memory = 
    }
    command {
        set -e 
        module load ${trimgaloreModule} && \
        if [ "${inputFastq2}x" == "x" ];then \
            trim_galore "${inputFastq1}" \
                --output_dir "$(basename $(basename ${outputFastq1} .fastq.gz) .fq.gz)"
                ln -sf $PWD/$(basename $(basename ${outputFastq1} .fastq.gz) .fq.gz)//reads_R1_val_1.fq.gz ${outputFastq1}
        else \
            trim_galore --paired "${inputFastq1}" "${inputFastq2}" \
              --output_dir "$(basename $(basename ${outputFastq1} .fastq.gz) .fq.gz)" 
            ln -sf $PWD/$(basename $(basename ${outputFastq1} .fastq.gz) .fq.gz)/reads_R1_val_1.fq.gz ${outputFastq1}
            ln -sf $PWD/$(basename $(basename ${outputFastq1} .fastq.gz) .fq.gz)//reads_R1_val_2.fq.gz ${outputFastq2}
        fi
    }

    output {
        File fastq1 = outputFastq1
        File? fastq2 = outputFastq2
    }

    runtime {
        memory: select_first([memoryGb * 1024,1024])

    }
}
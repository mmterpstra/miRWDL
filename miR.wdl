version 1.0

import "structs.wdl" as structs
import "tasks/common.wdl" as common
import "tasks/trimgalore.wdl" as trimgalore
import "tasks/fastqc.wdl" as fastqc

workflow fastqQuantWorkflow {
    input {
        #File fastq_input
        #Int mem_gb
	    #String? fastqc_mod = "FastQC/0.11.9-Java-11"
        #Int? fastqc_mem = 1
        #String? trimgalore_mod = "Trim_Galore/0.6.6-GCCcore-9.3.0-Python-3.8.2"
        #String? gatk_mod = "GATK/4.2.6.1-GCCcore-11.2.0-Java-11"
        File sampleJson
    }
    SampleConfig sampleConfig = read_json(sampleJson)
    
    scatter (sample in sampleConfig.samples) {
        #Array[ReadGroup] readgroups = sample.readgroups
        #Array[File] fastqList = flatten(readgroups["fastq1"])
        scatter (rg in sample.readgroups) {
            call common.CreateLink as getfastq1 {
            input:
                inputFile = rg.fastq1,
                outputPath = sample.name + "_" + rg.identifier + "_R1.fastq.gz"
            }
            call fastqc.FastQC as fastqc1 {
            input:
                inputFastq = getfastq1.link,
            }
            if (defined(rg.fastq2)) {
                call common.CreateLink as getfastq2 {
                input:
                    inputFile = select_first([rg.fastq2]),
                    outputPath = sample.name + "_" + rg.identifier + "_R2.fastq.gz"
                }
                call fastqc.FastQC as fastqc2 {
                input:
                inputFastq = getfastq2.link,
                }

            }
            #if (defined(rg.fastq2)) {
            call trimgalore.TrimGalore as adaptertrim {
                input:
                    inputFastq1 = getfastq1.link,
                    inputFastq2 = getfastq2.link,
                    outputFastq1 = getfastq1.link,
                    outputFastq2 = getfastq2.link,
                    memoryGb = 1
            }
            #}
            #if (!defined(rg.fastq2)) {
            #    call trimgalore.TrimGalore as adaptertrimPe {
            #        input:
            #            inputFastq1 = getfastq1.link,
            #            outputFastq1 = getfastq1.link,
            #            memoryGb = 1
            #    }
            #} 
        }       
        #call common.ConcatenateTextFiles as CombineSampleFastqGz{
        #    input:
        #        fileList=,
        #        combinedFilePath=sample.name + ".fastq.gz"
        #}

        call common.CollapseFastq as collapse {
            input:
                reads = select_first([adaptertrim.fastq1]),
                outputPrefix = sample.name,
                threeLetterName = sample.threeLetterName
        }
        
    }
}
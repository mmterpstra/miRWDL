version 1.0

import "structs.wdl" as structs
import "tasks/common.wdl" as common
import "tasks/trimgalore.wdl" as trimgalore
import "tasks/fastqc.wdl" as fastqc
import "tasks/mirdeep2.wdl" as mirdeep
import "tasks/multiqc.wdl" as multiqc


workflow fastqQuantWorkflow {
    input {
	    String fastqcModule = "FastQC/0.11.9-Java-11"
        String trimgaloreModule = "Trim_Galore/0.6.6-GCCcore-9.3.0-Python-3.8.2"
        String mirdeepModule = "miRDeep2/0.1.1-foss-2018b"
        String multiqcModule = "multiqc/1.12-GCCcore-11.3.0"
        File sampleJson
        File mirbaseHairpinFasta
        File mirbaseMatureFasta
    }

    SampleConfig sampleConfig = read_json(sampleJson)

    call mirdeep.ExtractMiRNAs as extractMiRNAs {
        input:
            mirdeepModule = mirdeepModule,
            inputMirbaseMatureFasta = mirbaseMatureFasta,
            inputMirbaseHairpinFasta = mirbaseHairpinFasta,
            speciesCode = "hsa,ebv,hcmv,jcv",
            extractedHairpinPrefix = "subset_hairpin",
            extractedMaturePrefix = "subset_mature"
    }
    scatter (sample in sampleConfig.samples) {
        #Array[ReadGroup] readgroups = sample.readgroups
        #Array[File] fastqList = flatten(readgroups["fastq1"])
        scatter (rg in sample.readgroups) {
            #linking for uniform filenames
            call common.CreateLink as getfastq1 {
            input:
                inputFile = rg.fastq1,
                outputPath = sample.name + "_" + rg.flowcell + "_" + rg.identifier + "_R1.fastq.gz"
            }
            call fastqc.FastQC as fastqc1 {
            input:
                inputFastq = getfastq1.link,
                fastqcModule = fastqcModule
            }
            if (defined(rg.fastq2)) {
                call common.CreateLink as getfastq2 {
                input:
                    inputFile = select_first([rg.fastq2]),
                    outputPath = sample.name + "_" + rg.flowcell + "_" + rg.identifier + "_R2.fastq.gz"
                }
                call fastqc.FastQC as fastqc2 {
                input:
                    inputFastq = getfastq2.link,
                    fastqcModule = fastqcModule
                }

            }

            #adapter trim these are assumend to be illumina adapters
            call trimgalore.TrimGalore as adaptertrim {
                input:
                    inputFastq1 = getfastq1.link,
                    inputFastq2 = getfastq2.link,
                    outputFastq1 = sample.name + "_" + rg.identifier + "_trim_R1.fastq.gz",
                    outputFastq2 = sample.name + "_" + rg.identifier + "_trim_R2.fastq.gz",
                    memoryGb = 1,
                    trimgaloreModule = trimgaloreModule
            }
            #}
            #if (!defined(rg.fastq2)) {
            #    call trimgalore.TrimGalore as adaptertrimSe {
            #        input:
            #            inputFastq1 = getfastq1.link,
            #            outputFastq1 = sample.name + "_" + rg.identifier + "_trim_R1.fastq.gz",
            #            memoryGb = 1
            #    }
            #} 
        }       
        #call common.ConcatenateTextFiles as CombineSampleFastqGz{
        #    input:
        #        fileList=,
        #        combinedFilePath=sample.name + ".fastq.gz"
        #}
        call fastqc.FastQCSample as fastqcSample1 {
            input:
                fastqcModule = fastqcModule,
                inputFastqGzs = getfastq1.link,
                outputPrefix = sample.name + "_R1",
        }
        if (defined(getfastq2.link)) {
            call fastqc.FastQCSample as fastqcSample2 {
            input:
                fastqcModule = fastqcModule,
                inputFastqGzs = select_all(getfastq2.link),
                outputPrefix = sample.name + "_R2",
            }
        }
        # Format conversion to collapsed mirdeep format for quantifier
        call common.CollapseFastq as collapse {
            input:
                reads = adaptertrim.fastq1,
                nextflex = true,
                outputPrefix = sample.name,
                threeLetterName = sample.threeLetterName
        }

        # Quantify against entire mirbase set
        call mirdeep.QuantifierSingleSample as quantify {
             input:
                mirdeepModule = mirdeepModule,
                inputCollapsedFasta = collapse.outputCollapsedFasta,
                inputMatureFasta = extractMiRNAs.outMatureFa,
                inputHairpinFasta = extractMiRNAs.outHairpinFa,
                outputPrefix = sample.name
        }
        
    }
    call mirdeep.MergeQuantifierOutputs as mergeOutputsQuantifier {
        input:
                inputLogs = quantify.outLog,
                inputTsvs = quantify.outTsv,
                inputCollapsedFasta = collapse.outputCollapsedFasta,
                outputPrefix = "quantifier_final"
    }
    #insert select_all(fastqcSample2.outZip),select_all(flatten(fastqc2.outZip)) when pe
    Array[File] files = flatten(flatten([fastqc1.outZip,adaptertrim.fastq1Log,[fastqcSample1.outZip]]))
    #Array[String?] optionalFiles = ["None"]
    #if (defined(select_first(fastqc2.outZip))) {
    #   Array[String] optionalFiles =  select_all(flatten(flatten([fastqc2.outZip,[fastqcSample2.outZip]])))
    #}
    #hmm passing variables as a task call will result in different resolution... gettin past the validation step.
    #Dont ask questions this thing below works for correct paired end data
    # extra flatten because option only accepts single 
    call multiqc.MultiQC as multiQC {
        input:
            multiqcModule = multiqcModule,
            files = files,
            optionalFiles = select_all(flatten(flatten([fastqc2.outZip,[fastqcSample2.outZip]])))
    }
}
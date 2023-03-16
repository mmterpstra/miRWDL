version 1.0

task QuantifierSingleSample {
    input {
        File inputCollapsedFasta
        File inputMirbaseHairpinFasta
        #File inputSampleConfig
        File inputMirbaseMatureFasta
        String outputPrefix = 'quantifier'
        Int? memoryGb = 1 + ceil(size(inputCollapsedFasta, "G"))*1
	    String mirdeepModule  = "mirdeep2/0.1.3-GCC-10.2.0-Perl-5.32.0"
        Int timeMinutes = 1 + ceil(size(inputCollapsedFasta, "G")) * 10
    }

    command {
        set -e
        module load ${mirdeepModule}
        
        quantifier.pl \
            -p "${inputMirbaseHairpinFasta}" \
            -m "${inputMirbaseMatureFasta}" \
            -r "${inputCollapsedFasta}" \
            -y "${outputPrefix}" \
            -P -W 2> "quantifier.stderr"
        
        cat quantifier.stderr > /dev/stderr
        
        #the stderr of quantifier is also the qc log for mapping results... so I catch the qc part by omitting the progress/warning messages
        grep -v 'Could not determine where \|Putting it to the 5p species hash\|creating pdf for' "quantifier.stderr" > \
            "quantifier_""$(basename ~{outputPrefix})"".log"
    }

    output {
        File outTsv = "miRNAs_expressed_all_samples_" + basename(outputPrefix, ".tsv") + ".csv"
        File outHtml = "expression_" + basename(outputPrefix, ".tsv") + ".html"
        File? outPdfs = "pdfs_" + basename(outputPrefix, ".tsv")
        File outLog = "quantifier_" + basename(outputPrefix, ".tsv") + ".log"
    }

    runtime {
        memory: select_first([memoryGb * 1024,1024])
        timeMinutes: timeMinutes
    }
}

task MergeQuantifierOutputs {
    input {
        Array[File] inputTsvs
        Array[File] inputLogs
        Array[File] inputCollapsedFasta
        String outputPrefix
        Int? memoryGb = 1
	    #String mirdeepModule  = "mirdeep2/0.1.3-GCC-10.2.0-Perl-5.32.0"
        Int timeMinutes = 30
    }

    command <<<
        set -e
        set -o pipefail
        
        touch ~{outputPrefix}.tsv
        for i in  ~{sep=" " inputTsvs}; do
            if [ ! -s ~{outputPrefix}.tsv ]; then 
                (echo -e "test\t\t";cut -f 1-3 $i )>> ~{outputPrefix}.tsv
            fi
            paste ~{outputPrefix}.tsv <(echo -e $(basename $i)"\t"; cut -f 5,6  $i )> ~{outputPrefix}.tmp.csv
            mv  ~{outputPrefix}.tmp.csv ~{outputPrefix}.tsv
        done

        (for i in  ~{sep=" " inputCollapsedFasta}; do 
            perl  -wne '
                BEGIN{our %lencounts;
                    our @h;
                    our @lens=(1..151)};
                    if($.%2==1){ 
                        @h = split("_x",$_);
                    }else{
                        chomp;
                        $lencounts{length($_)}+=$h[1];
                    }
                    END{
                        #print $ARGV."\n";
                        print join("\t",("File",@lens))."\n";
                        my @counts;
                        for my $key (@lens){
                            if(defined($lencounts{$key})){
                                push(@counts,$lencounts{$key});
                            }else{
                                push(@counts,0)
                            }
                        };
                        print join("\t",($ARGV,@counts))."\n";
                    }' $i 
        done) > ~{outputPrefix}".collapsedcounts.log"

        (

            echo "### QC metrics of miRNA quantification workflow"
            cat << END >/dev/stdout| perl -wpe 's/^ +//'
            
            Per sample metrics are present in this file. Other metrics are 
            present in the in the 'multiqc.html' file.
            
            The core workflow looks like:
            
            adapter trimming (TrimGalore) -> collapse (scripting) -> 
            quantification (quantifier.pl from mirdeep2). 
            
            The adapter trimming is performed on a fastq input level and those 
            are merged during collapse of data. Depending on the type of data 
            N4 barcodes are clipped during collapse.

            FastQC and MultiQC are also used for qualitity control.

            > MultiQC
            > MultiQC: Summarize analysis results for multiple tools and samples in a single report
            > Philip Ewels, Måns Magnusson, Sverker Lundin and Max Käller
            > Bioinformatics (2016)
            > doi: 10.1093/bioinformatics/btw354
            > PMID: 27312411
         
            > Cutadapt (used in TrimGalore)
            > MARTIN, Marcel. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, [S.l.], v. 17, n. 1, p. pp. 10-12, may 2011. ISSN 2226-6089. Available at: <https://journal.embnet.org/index.php/embnetjournal/article/view/200>. Date accessed: 14 feb. 2023. doi:https://doi.org/10.14806/ej.17.1.200. 

            > miRDeep2
            > Friedländer MR, Mackowiak SD, Li N, Chen W, Rajewsky N. miRDeep2 accurately identifies known and hundreds of novel microRNA genes in seven animal clades. Nucleic Acids Res. 2012 Jan;40(1):37-52. doi: 10.1093/nar/gkr688. Epub 2011 Sep 12. PMID: 21911355; PMCID: PMC3245920.

END
            echo ""
            echo "### Read length counts"
            echo ""
            echo "Per sample counts of an observed read length."
            
            perl -wpe 'chomp;
                        $_="| ".$_." |\n";
                        s/\t/ | /g;
                        if($. == 1){
                            print $_ ;
                            s/\| .*\/\| /|/g;
                            s/\.md\.fa//g;
                            s/[\w\\\/\.]+/---/g
                        }' ~{outputPrefix}".collapsedcounts.log" | column -t -s\ 
        )>  ~{outputPrefix}".md"
        
        (
            echo ""
            echo "### Quantifier mapping yield"
            echo ""
            echo "Per sample metrics of reads mapped against mirbase."
            echo ""

            (   for i in ~{sep=" " inputLogs}; do
                    echo ""
                    head -n 20 $i | \
                        grep 'Expressed miRNAs are written to expression_analyses/expression_analyses' | \
                        perl -wpe 's!Expressed miRNAs are written to expression_analyses\/expression_analyses_|/miRNA_expressed.csv!!g'
                    echo ""
                    head -n 20  $i | \
                        grep '#desc' | \
                        perl -wpe 'chomp; s/^#|$|\t/ | /g;$_.="\n";s/^ +//g;print $_; s/[a-z%]/-/g;'
                    head -n 20  $i | \
                        grep '#desc' -A 2 | \
                        tail -n 1 | \
                        perl -wpe 'chomp; s/^|$|\t|:/ | /g;s/ +/ /g;s/^ +//g;$_.="\n";'
                done
            ) | column -t -s\ 
        )>>  ~{outputPrefix}".md"

    >>>
    output {
        File outTsv = basename(outputPrefix, ".tsv") + ".tsv"
        File outLog = basename(outputPrefix, ".tsv") + ".md"
        File outCollapseLog = basename(outputPrefix, ".tsv") + ".collapsedcounts.log"
    }
    runtime {
        memory: select_first([memoryGb * 1024,1024])
        timeMinutes: timeMinutes
    }
}

task Quantifier {
    input {
        Array[File] inputCollapsedFasta
        File inputMirbaseHairpinFasta
        File inputSampleConfig
        File inputMirbaseMatureFasta
        String speciesCode = 'hsa'
        String outputPrefix = 'quant'
        Int? memoryGb = 1 
        #+ ceil(size(inputCollapsedFasta, "G"))*1
	    String mirdeepModule  = "mirdeep2/0.1.3-GCC-10.2.0-Perl-5.32.0"
        Int timeMinutes = 1
        # + ceil(size(inputCollapsedFasta, "G")) * 10
    }

    command <<<
        set -e
        module load ~{mirdeepModule}
        cat ~{sep=" " inputCollapsedFasta} > "~{outputPrefix}.all_reads.md_collapse.fa"
        #this probably works fine for the first 2600 or so samples 
        paste  <(ls ~{sep=" " inputCollapsedFasta}) \
            <(ls ~{sep=" " inputCollapsedFasta} | \
                python3 -c "import sys;import math;[sys.stdout.write(chr(65+math.floor((count+1)/99)) + ':02d'.format((count+1)%99)) for count,line in enumerate(sys.stdin)]") \
            > "~{outputPrefix}.sample.config"
        
        quantifier.pl -c "~{inputSampleConfig}" \
            -p "~{inputMirbaseHairpinFasta}" \
            -m "~{inputMirbaseMatureFasta}" \
            -r "~{inputCollapsedFasta}" \
            -y "~{outputPrefix}" \
            -t "~{speciesCode}" \
            -P -W 2> "quantifier.stderr"

        cat quantifier.stderr > /dev/stderr
        #the stderr of quantifier is also the qc log for mapping results... so I catch the qc part by omitting the progress/warning messages
        grep -v 'Could not determine where \|Putting it to the 5p species hash\|creating pdf f or' "quantifier.stderr" > \
            "quantifier_""~{outputPrefix}"".log"
    >>>

    output {
        #String fastqcDir =  basename(inputCollapsedFasta, ".fastq.gz")
        #ile out_html1 = basename(inputCollapsedFasta, ".fastq.gz") + "_fastqc.html"
        File outCsv = "miRNAs_expressed_all_samples_" + basename(outputPrefix, ".tsv") + ".csv"
        File outHtml = "expression_" + basename(outputPrefix, ".tsv") + ".html"
        File? outPdfs = "pdfs_" + basename(outputPrefix, ".tsv")
        File outLog = "quantifier_"+ basename(outputPrefix, ".tsv") + ".log"
    }

    runtime {
        memory: select_first([memoryGb * 1024,1024])
        timeMinutes: timeMinutes
    }
}

#task to download mirbase?

task ExtractMiRNAs {
    input {
        File inputMirbaseHairpinFasta
        File inputMirbaseMatureFasta
        String extractedHairpinPrefix
        String extractedMaturePrefix
        String speciesCode = 'hsa'
        Int? memoryGb = 1 
        #+ ceil(size(inputCollapsedFasta, "G"))*1
	    String mirdeepModule  = "mirdeep2/0.1.3-GCC-10.2.0-Perl-5.32.0"
        Int timeMinutes = 15
        # + ceil(size(inputCollapsedFasta, "G")) * 10
    }

    command <<<
        set -e
        module load ~{mirdeepModule}
        
        extract_miRNAs.pl "~{inputMirbaseHairpinFasta}" "~{speciesCode}" > "~{extractedHairpinPrefix}.fa"
        extract_miRNAs.pl "~{inputMirbaseMatureFasta}" "~{speciesCode}" mature > "~{extractedMaturePrefix}.fa"
        

        
    >>>

    output {
        File outHairpinFa = extractedHairpinPrefix + ".fa"
        File outMatureFa = extractedMaturePrefix + ".fa"
    }

    runtime {
        memory: select_first([memoryGb * 1024,1024])
        timeMinutes: timeMinutes
    }
}

version 1.0

task QuantifierSingleSample {
    input {
        File inputCollapsedFasta
        File inputHairpinFasta
        #File inputSampleConfig
        File inputMatureFasta
        String outputPrefix = 'quantifier'
        Int? memoryGb = 2 + ceil(size(inputCollapsedFasta, "G"))*1
	    String mirdeepModule
        #the files below rarely exceed 1-2gb 
        Int timeMinutes = 1 + ceil(size(inputHairpinFasta, "G")) * 30 + ceil(size(inputCollapsedFasta, "G")) * 90
    }

    command {
        set -e
        
        module --ignore-cache load ${mirdeepModule}
        mkdir -p "expression_analyses/expression_analyses_${outputPrefix}_test/"
        #this command below breaks on execution nodes idk why
        (cd "expression_analyses/expression_analyses_${outputPrefix}_test/" && python3 $EBROOTBOWTIE/bin/bowtie-build ${inputHairpinFasta} ./miRNA_precursor)
        (ldd   $(which bowtie-build-l); ldd   $(which bowtie-build-s))
        module --ignore-cache load ${mirdeepModule}
        quantifier.pl \
            -p "${inputHairpinFasta}" \
            -m "${inputMatureFasta}" \
            -r "${inputCollapsedFasta}" \
            -y "${outputPrefix}" \
            -P -W 2> "quantifier.stderr"
        
        if [ -e "pdfs_${outputPrefix}" ]; then 
            zip -ru "pdfs_${outputPrefix}.zip" "pdfs_${outputPrefix}"
        fi

        cat quantifier.stderr > /dev/stderr
        
        #the stderr of quantifier is also the qc log for mapping results... so I catch the qc part by omitting the progress/warning messages
        grep -v 'Could not determine where \|Putting it to the 5p species hash\|creating pdf for' "quantifier.stderr" > \
            "quantifier_""$(basename ~{outputPrefix})"".log"
    }

    output {
        File outTsv = "miRNAs_expressed_all_samples_" + basename(outputPrefix, ".tsv") + ".csv"
        File outHtml = "expression_" + basename(outputPrefix, ".tsv") + ".html"
        File? outPdfs = "pdfs_" + basename(outputPrefix, ".tsv") + ".zip"
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
        Array[File] inputFastqcZips
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
            SAMPLENAME="$(echo $i | perl -wpe 's/^.*\/miRNAs_expressed_all_samples_|\.csv$//g')"
            if [ ! -s ~{outputPrefix}.tsv ]; then 
                (cut -f 1-3 $i )>> ~{outputPrefix}.tsv
            fi
            paste ~{outputPrefix}.tsv <(echo -e $SAMPLENAME"\t"$SAMPLENAME"(norm)"; cut -f 5,6  $i | tail -n+2)> ~{outputPrefix}.tmp.csv
            mv  ~{outputPrefix}.tmp.csv ~{outputPrefix}.tsv
        done

        (for i in  ~{sep=" " inputCollapsedFasta}; do 
            #note on @h = split("_x",$_):
            #mirdeep collapsed fasta format header assumed
            #(>\w[\w\d]{2})_(\d+)_x(\d+) $1: 3 letter sample descriptor. $2: collapsed readid. $3: collapsed read count.
            perl  -wne '
                BEGIN{our %lencounts;
                    our @h;
                    our @lens=("Total",1..151)};
                if($.%2==1){ 
                    @h = split("_x",$_);
                }else{
                    chomp;
                    my $count = 1;
                    $count = $h[1] if(defined($h[1]));
                    if(defined($lencounts{length($_)})){
                        $lencounts{length($_)}+=$count;
                    }else{
                        $lencounts{length($_)}=$count;
                    }
                    if(defined($lencounts{"Total"})){
                        $lencounts{"Total"}+=$count;
                    }else{
                        $lencounts{"Total"}=$count;
                    }
                }
                END{
                    #print $ARGV."\n";
                    print join(",",("File",@lens))."\n";
                    my @counts;
                    for my $key (@lens){
                        if(defined($lencounts{$key})){
                            push(@counts,$lencounts{$key});
                        }else{
                            push(@counts,0)
                        }
                    };
                    chomp($ARGV);
                    print join(",",($ARGV,@counts))."\n";
                }' $i 
        done) | perl -wpe 's/^.*\/call-mergeOutputsQuantifier\/inputs\/[-\d]+\///g' > ~{outputPrefix}".collapsedcounts.log"

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
            echo "### Read counts input data"
            echo ""
            echo "For each sample the input reads."
            echo ""
            (echo -e "Sample\tReads";
            for i in  ~{sep=" " inputFastqcZips}; do 
                echo -ne "$(basename "$i" '_R1_fastqc.zip' )""\t"
                unzip -p "$i"  \*/fastqc_data.txt | head -n 10 | perl -wne 'print if s/Total Sequences\s+//'
            done)| \
                perl -wpe 'chomp;
                    s/^|$|\t|:/ | /g;
                    s/ +/ /g;
                    s/^ +| +$//g;
                    if($. == 1){
                        print $_."\n";
                        s/\| .*\/\| /|/g;
                        s/\.md\.fa//g;
                        s/[^\s|]+/---/g;
                    }
                    $_.="\n";' | \
                column -t -s\   
            
            echo ""

            echo ""
            echo "### Read length counts post trim"
            echo ""
            echo "Readcounts for each read length after trimming for each sample. These should show a peak around miRNA length (22 bases)."
            
            perl -wne 'chomp;
                        $_="| ".$_." |\n";
                        s/,/ | /g;
                        if($. == 1){
                            print $_ ;
                            s/\| .*\/\| /|/g;
                            s/\.md\.fa//g;
                            s/[\w\\\/\.]+/---/g;
                            print $_

                        }else{
                            if( not(m/\| File\s*\|/) ){
                                print $_;
                            }
                        }' ~{outputPrefix}".collapsedcounts.log" | column -t -s\ 
        )>  ~{outputPrefix}".md"
        
        (
            echo ""
            echo "### Quantifier mapping yield"
            echo ""
            echo "Per sample metrics of reads mapped against mirbase. Total reads is reads remaining after adapter trim and after import which discards N containing sequences"
            echo ""

            (   for i in ~{sep=" " inputLogs}; do
                    echo ""
                    #samplename
                    SAMPLENAME=$(head -n 20 $i | \
                        grep 'Expressed miRNAs are written to expression_analyses/expression_analyses' | \
                        perl -wpe 's!Expressed miRNAs are written to expression_analyses\/expression_analyses_|/miRNA_expressed.csv!!g;chomp')
                    echo ""
                    if [ ! -e "~{outputPrefix}.md.header" ] ;then
                        #header and (default) markdown table alignment
                        head -n 20  $i | \
                            grep '#desc' | \
                            perl -wpe 'chomp; s/^#|$|\t/ | /g;$_.="\n";s/^ +//g;print $_; s/[a-z%]/-/g;'
                        touch "~{outputPrefix}.md.header"
                    fi

                    #sample info
                    (echo -en "$SAMPLENAME\t"
                        head -n 20  $i | \
                        grep '#desc' -A 2 | \
                        tail -n 1 | \
                        cut -f 2- -d\: ) | \
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
        module --ignore-cache load ~{mirdeepModule}
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
	    String mirdeepModule
        Int timeMinutes = 15
        # + ceil(size(inputCollapsedFasta, "G")) * 10
    }

    command <<<
        set -e
        module --ignore-cache load ~{mirdeepModule}
        
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
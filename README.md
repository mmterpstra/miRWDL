# miRWDL
miRWDL scalable mirna quantification on the wdl with a default configuration for slurm

### Example
Example how to run at test/test.sh

### How to install

Installation requirements 

 - Fastqc
 - TrimGalore
    - Cutadapt
 - miRDeep2
    - Bowtie(1)

Installation by easybuild (preferred or by configuring docker images)

Tested versions
   - FastQC FastQC/0.11.9-Java-11 
   - TrimGalore Trim_Galore/0.6.6-GCCcore-9.3.0-Python-3.8.2
   - miRDeep2 mirdeep2/0.1.3-GCC-10.2.0-Perl-5.32.0

### Future work

Quantify other small RNA subgroups taking the groups by exclusion of ensembl genetypes.

```
cat  /path/to/ftp.ensembl.org/pub/release-99/gtf/Homo_sapiens.GRCh38.99.gtf | grep -v 'unprocessed_pseudogene\|lncRNA\|"unprocessed_pseudogene"\|miRNA\|"transcribed_processed_pseudogene"\|"processed_pseudogene"\|"protein_coding"\|"unitary_pseudogene"\|"transcribed_unitary_pseudogene"\|"polymorphic_pseudogene"\|"translated_processed_pseudogene"\|"IG_V_pseudogene"\|"IG_C_gene"\|"IG_J_gene"\|"IG_V_gene"\|"IG_C_pseudogene"\|"IG_D_gene"\|"IG_J_pseudogene"\|"IG_pseudogene"\|"TR_[CDJV]_gene"\|"TR_[JV]_pseudogene"\|"pseudogene"' | perl -wpe 's/.*gene_biotype "(.*?)";.*/$1/g;' | sort -u 
#!genebuild-last-updated 2019-08
#!genome-build-accession NCBI:GCA_000001405.28
#!genome-build GRCh38.p13
#!genome-date 2013-12
#!genome-version GRCh38
misc_RNA
Mt_rRNA
Mt_tRNA
ribozyme
rRNA
rRNA_pseudogene
scaRNA
scRNA
snoRNA
snRNA
sRNA
TEC
vaultRNA
#amount of genes -5
cat /path/to/ftp.ensembl.org/pub/release-99/gtf/Homo_sapiens.GRCh38.99.gtf | grep -v 'unprocessed_pseudogene\|lncRNA\|"unprocessed_pseudogene"\|miRNA\|"transcribed_processed_pseudogene"\|"processed_pseudogene"\|"protein_coding"\|"unitary_pseudogene"\|"transcribed_unitary_pseudogene"\|"polymorphic_pseudogene"\|"translated_processed_pseudogene"\|"IG_V_pseudogene"\|"IG_C_gene"\|"IG_J_gene"\|"IG_V_gene"\|"IG_C_pseudogene"\|"IG_D_gene"\|"IG_J_pseudogene"\|"IG_pseudogene"\|"TR_[CDJV]_gene"\|"TR_[JV]_pseudogene"\|"pseudogene"' | cut -f 1,4,5,7 | sort -u 
6823
#convert into sequence fasta for quantification
#first create bed
grep -v '#' /path/to/ftp.ensembl.org/pub/release-99/gtf/Homo_sapiens.GRCh38.99.gtf | grep -v 'unprocessed_pseudogene\|lncRNA\|"unprocessed_pseudogene"\|miRNA\|"transcribed_processed_pseudogene"\|"processed_pseudogene"\|"protein_coding"\|"unitary_pseudogene"\|"transcribed_unitary_pseudogene"\|"polymorphic_pseudogene"\|"translated_processed_pseudogene"\|"IG_V_pseudogene"\|"IG_C_gene"\|"IG_J_gene"\|"IG_V_gene"\|"IG_C_pseudogene"\|"IG_D_gene"\|"IG_J_pseudogene"\|"IG_pseudogene"\|"TR_[CDJV]_gene"\|"TR_[JV]_pseudogene"\|"pseudogene"' | cut -f 1,4,5,7,9- |  sort -k1,1g -k2,2n -k3,3n -k4g -s$(echo -e "\t") -u | perl -F'/\t/' -wane '$F[1] +=-1;$F[4] =~ s/.*"(ENS[GT]\d+?)".*gene_name "(\S*?)";.*/${1}_${2}/g;chomp $F[4];print join("\t",(@F[0..2],$F[4],"0",$F[3]))."\n"' | head
MT	576	647	ENSG00000210049_MT-TF	0	+
MT	647	1601	ENSG00000211459_MT-RNR1	0	+
MT	1601	1670	ENSG00000210077_MT-TV	0	+
MT	1670	3229	ENSG00000210082_MT-RNR2	0	+
MT	3229	3304	ENSG00000209082_MT-TL1	0	+
MT	4262	4331	ENSG00000210100_MT-TI	0	+
MT	4328	4400	ENSG00000210107_MT-TQ	0	-
MT	4401	4469	ENSG00000210112_MT-TM	0	+
MT	5511	5579	ENSG00000210117_MT-TW	0	+
MT	5586	5655	ENSG00000210127_MT-TA	0	-
#convert into unique sequence noncoding rna fasta using bedtools
grep --color=auto -v '#' /path/to/ftp.ensembl.org/pub/release-99/gtf/Homo_sapiens.GRCh38.99.gtf | \
grep --color=auto -v 'unprocessed_pseudogene\|lncRNA\|"unprocessed_pseudogene"\|miRNA\|"transcribed_processed_pseudogene"\|"processed_pseudogene"\|"protein_coding"\|"unitary_pseudogene"\|"transcribed_unitary_pseudogene"\|"polymorphic_pseudogene"\|"translated_processed_pseudogene"\|"IG_V_pseudogene"\|"IG_C_gene"\|"IG_J_gene"\|"IG_V_gene"\|"IG_C_pseudogene"\|"IG_D_gene"\|"IG_J_pseudogene"\|"IG_pseudogene"\|"TR_[CDJV]_gene"\|"TR_[JV]_pseudogene"\|"pseudogene"' | \
cut -f 1,4,5,7,9- | \
sort -k1,1g -k2,2n -k3,3n -k4g -s$(echo -e "\t") -u | \
perl -F'/\t/' -wane '$F[1] +=-1;$F[4] =~ s/.*"(ENS[GT]\d+?)".*gene_name "(\S*?)";.*/${1}_${2}/g;$F[0] = "chr".$F[0];$F[0] =~ s/chr([KG][LI]\d+)\.(\d)/chrUn_${1}v${2}/g if($F[0] =~ m/chrGL\d+|chrKI\d+/);$F[0]="chrM" if($F[0] eq "chrMT");chomp $F[4];print join("\t",(@F[0..2],$F[4],"0",$F[3]))."\n"' | \
perl -wne 'chomp; print $_."\t".`ml BEDTools; echo "$_"| bedtools getfasta -fi /path/to/ftp.broadinstitute.org/bundle/bundle17jan2020/hg38/Homo_sapiens_assembly38.fasta -bed /dev/stdin | \
tail -n1 | perl -wpe 'chomp'`."\n";' | \
perl -wane 'print ">".$F[3]." ".$F[0]."-".$F[1].":".$F[2]."\n".$F[6]."\n" if($F[6] ne "");' | \
perl  -wpe 'if($.%2 == 1){chomp;$_ .= "\t";}' | \
perl -F"\t" -wane 'use Data::Dumper; BEGIN{our $fastaToHeaderList;%{$fastaToHeaderList} }; if(defined($fastaToHeaderList -> {$F[1]})){push(@{$fastaToHeaderList -> {$F[1]}},substr($F[0],1))}else{$fastaToHeaderList -> {$F[1]} = [substr($F[0],1)]}; warn Dumper($fastaToHeaderList) if $. == 100 ;END{ for my $fa (keys(%{$fastaToHeaderList})){print ">".join(" ", @{$fastaToHeaderList -> {$fa}})."\n"; print $fa}};'> /path/to/ftp.ensembl.org/pub/release-99/gtf/ncRNA_Homo_sapiens.GRCh38.99.fasta
```